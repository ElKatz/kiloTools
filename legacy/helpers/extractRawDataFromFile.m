function [raw] = extractRawDataFromFile(fullPathData, dataFileType, chString, opts)
%   [raw] = extractRawDataFromFile(fullPathData, dataFileType, chString, opts)
%

switch dataFileType
    
    case 'plx'
        %% create a list of all ad continuous channel names in cell array:
        [nCh, chName]     = plx_adchan_names(fullPathData);
        chNameList = cell(nCh,1);
        for ii = 1:nCh
            chNameList{ii} = chName(ii,:);
        end
        
        %% get data for channel defined by 'channelString':
        % get indices for the spike channels & lfp channels:
        idxCh = getIdx_forChannelType(chNameList, chString);
        
        % get number of spikes counts per ad channel & get indices for
        % those that have data:
        [~, samplecounts] = plx_adchan_samplecounts(fullPathData);
        idxHaveData = samplecounts~=0;
        
        % good indices are those that are spk channels & have data:
        idxGood   = idxCh & idxHaveData;
        
        % good spk channles names:
        chNameGood   = chNameList(idxGood);
        
        %% build data matrix 'samples' of size [nChannels, nSamples]:
        
        % nChannels & nSamples:
        nChannels   = numel(chNameGood);
        tmp         = samplecounts(idxGood);
        nSamples    = tmp(1); % taking the number of samples in first spk channel. Rest are identical.
        
        % prealloc:
        samples     = zeros(nChannels, nSamples, 'int16');
        tChRead     = nan(nChannels,1); % time keeping
        
        % gotta map out indices to plxeon's ad channel numbers:
        [~,   chNumber]   = plx_ad_chanmap(fullPathData);
        chNumberGood = chNumber(idxGood);
        
        % get it:
        fprintf('\tGetting data from %d %s channels...', nChannels, chString)
        for iCh = 1:nChannels
            [adfreq, n, ts, fn, ad] = plx_ad(fullPathData, chNumberGood(iCh)); % retunrs signal in miliVolts
            if n>0
                % data matrix 'samples':
                samples(iCh,1:n) = int16(ad);
            else
                disp('wtf?!')
                keyboard
            end
            % for first channel, save the ts & fn:
            if iCh==1
                ts1 = ts;
                fn1 = fn;
            end
        end
        
        %% get timing of ad data:
        % only doing this for the last channel because acquisition
        % time stamps & recording fragments are equal accross channels.
        fptinf('\tBuilding time vector. This step takes some time. My apologies')
        tSamples = [];
        for ii = 1:numel(ts1)
            tt = linspace(ts1(ii), ts1(ii)+fn1(ii)./adfreq, fn1(ii));
            tSamples = [tSamples(:); tt(:)];
        end
        fprintf('\tdone!')
        % GOTTA OPTIMIZE! linspace so slow....
        
        %% get strobed event time stamps and strobe values
        strobeChNumber = 257; % for amar rig, it is 257. might not be the case on other rigs.
        [~, tEvent, vEvent] = plx_event_ts(fullPathData, strobeChNumber);
        

        
    case 'matFromMpx'
        
        % * might need to include a step where I identify whether the mat
        % was converted in struct form or not.
        mpx     = load(fullPathData);
        chNameList  = {mpx.Channel_ID_Name_Map.Name};
        
        %% get data for channel defined by 'channelString':
        % get indices for the spike channels & lfp channels:
        idxCh = getIdx_forChannelType(chNameList, chString);
        idxGood = idxCh; % should I check for anything else? is there an 'idxGood' to make? check if all are zero otr someting? apply theo zero xhexk to plx too?
        
        % good spk channles names:
        chNameGood   = chNameList(idxGood);
        
        % nChannels:
        nChannels   = numel(chNameGood);
        
        % enforce equal sampling rate accross channels:
        mpx = enforceEqualSmplingRateAcrossChannels(mpx, chString);
        
        % Now that they are all the same Fs, taking the number of
        % samples in first spk channel:
        nSamples = numel(mpx.(chNameGood{1}).Samples);
        
        % build data matrix 'samples' of size [nChannels, nSamples]:
        samples     = zeros(nChannels, nSamples, 'int16');
        
        fprintf('\n\tGetting data from %d %s channels...', nChannels, chString)
        for iCh = 1:nChannels
            samples(iCh,:) = mpx.(chNameGood{iCh}).Samples;
        end
        
        % if you wish to only segemnt good epochs from the data (i.e. times
        % where animal was doing the task), by all means:
        if isfield(opts, 'performSegmentation') && opts.performSegmentation
            [samples, tSamples] = segmentSamples_mpx(samples, mpx, chString);
        else
            tSamples = linspace(mpx.(chNameGood{iCh}).TimeBegin, mpx.(chNameGood{iCh}).TimeEnd, numel(mpx.(chNameGood{iCh}).Samples));
        end
        
        % get event timestmpas (ts) and strobe values (sv):
        if isfield(mpx, 'CInPort_001')

            if isfield(mpx, 'SF_Sampling_Rate')
                tEvent  = mpx.CInPort_001.Samples(1,:) ./ (mpx.SF_Sampling_Rate*1e3);
            elseif isfield(mpx, 'SF_KHz')
                tEvent  = mpx.CInPort_001.Samples(1,:) ./ (mpx.SF_KHz*1e3);
            else
                error('couldn''t find a smapling right for the alphaLab Stream Format (SF) cause alphaLab suck')
            end
            vEvent  = mpx.CInPort_001.Samples(2,:);
        else
            tEvent = [];
            vEvent = [];
        end
        
end

%% output struct:
raw.samples     = samples;
raw.tSamples    = tSamples;
raw.tEvent      = tEvent;
raw.vEvent      = vEvent;
raw.chNameList  = chNameGood;

%% visualize
% 
% figure, 
% plot(tSamples, samples(1,:));
% hold on
% plot(tEvent, ones(1,numel(tEvent)), 'o')
% 
