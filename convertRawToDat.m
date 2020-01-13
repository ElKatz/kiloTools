function [samples, datPath, combo] = convertRawToDat(rawFullPath, opts)
%   [samples, datPath] = convertRawToDat(rawFullPath, opts)
%
% Converts continuous data from raw ephys file to matrix 'samples' of size:
% [nSamples, nChannels] and saves as binary .dat file in same folder as
% ephys file (unless specified otherwise in opts)
% INPUT:
%   rawFullPath  - Optional. full path to raw ephys data file.
%                   This file can be 'plx', 'mpx', 'oe', whatever, as long
%                   as it contains the continuous voltage traces.
%
% !!!!! AT THE MOMENT, THIS FUNCTION ONLY DEALS WITH plx/pl2 FILES !!!!!
%
%   opts - (optional) struct of options:
%   .outputFolder - full path to folder to save dat file (default: same
%                   folder as raw file)
%   .commonAverageReferencing - will subtract average over channels
%   .removeArtifacts - enter at own risk!
%   .specificChannels - user can select which plexon channels to use for
%                       conversion. remember, this must be in
%                       plexon-numbering, eg SPKC1 is usually ch num 65.
%   .plotProbeVoltage - if true, spits out a figure of the probe's voltage
%
% OUTPUT:
%   samples - [nChannels, nSamples] consisting of all continuous data
%   datPath - path to the .dat file


% 2do:
% generalize the identificaiton of continuous channel strings so that it
% works on plexon, alphaLab, openEphys etc...
%
% generalize to multiple file formats. vet.


dbstop if error

%% paths:
addPathsForSpikeSorting;

%% data file & folder names:

if ~exist('rawFullPath', 'var')
    [rawFileName, rawFolder] = uigetfile('*.*', 'Select files for conversion', '~/Dropbox/Code/spike_sorting/');
else
    [rawFolder, rawFileName, rawFileType]  = fileparts(rawFullPath);
    rawFileName = [rawFileName rawFileType];
    rawFileType = rawFileType(2:end);
end

% full path to plx file:
rawFullPath = fullfile(rawFolder, rawFileName);

% datasetname:
dsn = rawFileName(1:end-4);

%% use optional arguements and/or set defaults:
% init:
if ~exist('opts', 'var')
    opts = struct;
end

% output folder:
if ~isfield(opts, 'outputFolder')
    opts.outputFolder = fullfile(rawFolder, 'kiloSorted');
    if ~exist(opts.outputFolder, 'dir')
        mkdir(opts.outputFolder);
    end
end

%% options:
if ~isfield(opts, 'commonAverageReferencing')
    opts.commonAverageReferencing = false;
end

% remove artifacts
if ~isfield(opts, 'removeArtifacts')
    opts.removeArtifacts = false;
end

%% file names for .dat file (EPHYS) & .mat file (Timestamps and info):

% EPHYS: dat file named after dsn:
datPath = fullfile(opts.outputFolder, [dsn '.dat']);
% if a .dat file already exists delete it so that new file is so fresh and
% so clean clean
if exist(datPath, 'file')
    delete(datPath)
end

%% begin conversion:

if ~exist('rawFileType', 'var')
        rawFileType = 'pl2';
end
    
disp('--------------------------------------------------------------')
fprintf('Performing conversion of %s\n', dsn)
disp('--------------------------------------------------------------')
tic
% Different file types require different code to extract goodies. Each
% filetype (e.g. plx, mpx, etc.) gets its own case in this switch loop:
tStart = tic;
switch rawFileType
    
    case {'plx', 'pl2'}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% pl plx plx plx plx plx plx plx plx plx plx plx plx plx %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % create a list of all ad continuous channel names in cell array:
        [nCh, adChName]     = plx_adchan_names(rawFullPath);
        chNameList = cell(nCh,1);
        for ii = 1:nCh
            chNameList{ii} = adChName(ii,:);
        end
        idxSpkCh      = false(numel(chNameList),1);
        
        % if user provided specific channels to use for conversion, take
        % them:
        if isfield(opts, 'specificChannels') && ~isempty(opts.specificChannels)
            idxSpkCh(opts.specificChannels) = true;
        else
            % otherwise, figure out which spike-continuous channels have
            % data and grab'em:
            
            % get indices for the spike channels contunous ("SPKC")
            spkChannelStr = 'SPKC';
            spkChannelStr2 = 'CSPK'; % for alphaOmega converted mpx to plx files....
            for iCh = 1:numel(chNameList)
                if ~isempty(strfind(chNameList{iCh}, spkChannelStr)) || ~isempty(strfind(chNameList{iCh}, spkChannelStr2))
                    idxSpkCh(iCh) = true;
                else
                    idxSpkCh(iCh) = false;
                end
            end
        end
        
        % get number of spikes counts per ad channel and get use only 
        % those that have data:
        [~, samplecounts] = plx_adchan_samplecounts(rawFullPath);
        idxDataCh = samplecounts~=0;
        % get indices for channels that are both spk channels & have data:
        idxGoodCh =  idxSpkCh & idxDataCh;
        
        % nChannels & nSamples:
        nChannels   = sum(idxGoodCh);
        tmp         = samplecounts(idxGoodCh);
        nSamples    = tmp(1); % taking the number of samples in first spk channel. Rest are identical.
        
        % build data matrix 'samples' of size [nChannels, nSamples]:
        samples     = zeros(nChannels, nSamples, 'int16');
        tChRead     = nan(nChannels,1); % time keeping
        % gotta map out indices to plxeon's ad channel numbers:
        [~,   adChNumber]   = plx_ad_chanmap(rawFullPath);
        spkChNumber = adChNumber(idxGoodCh);
        fprintf('%0.1fs: Getting data from %0.0d spike channels!\n', toc, sum(idxGoodCh))
        hWait = waitbar(0, 'Converting channels...');
        for iCh = 1:nChannels
            tChRead(iCh) = toc;
            fprintf('\t%0.1fs: read channel #%0.0d \n', tChRead(iCh), spkChNumber(iCh));
            % data matrix 'samples':
            samples(iCh,:) = adReadWrapper(rawFullPath, spkChNumber(iCh)); % returns signal in miliVolts
            waitbar(iCh/nChannels, hWait, ['Converting channel ' num2str(iCh) ' of ' num2str(nChannels)]);
        end
        close(hWait)
        
        if isfield(opts, 'plotProbeVoltage') && opts.plotProbeVoltage
            plot_probeVoltage(samples, 1e3);
        end
        %%  extract timing information from raw file in "real" time

        %
        % 'tsMap' has a timestamp for every sample recorded. This will be a
        % vector of size nSamples. tsMap is used to convert from the spike 
        % index output of kiloSort (these are simply integers that indicate 
        % which sample number each spike occurred at), to time (in seconds)
        % relative to the beginning of the ephys recording.
        % This is needed because the event time stamps (evTs) from the raw
        % file are in same relative time (also in seconds).
        
        % get timestamps start values (tsStartVals) at start of each fragment:
        disp('Getting plexon timestamps for ad samples');
        
        % must read in a spike channel to construct the "timestamp map" from
        % samples (kilosort) to time in seconds.
        switch rawFileType
            case 'plx'
                [ad.ADFreq, ~,  ad.FragTs, ad.FragCounts] = plx_ad(rawFullPath, 'SPKC01');
            case 'pl2'
                ad = PL2Ad(rawFullPath, 'SPKC01');
        end
        % place to store the "map" from samples to seconds.
        sampsToSecsMap = zeros(sum(ad.FragCounts),1);

        % sample duration
        sampDur = 1/ad.ADFreq;

        % how many fragments of recording?
        nFrags = length(ad.FragTs);
        currentSample = 1;
        for i = 1:nFrags
            chunkIndex = currentSample:(currentSample + ad.FragCounts(i) - 1);
            timeStamps = ad.FragTs(i) + (0:(ad.FragCounts(i)-1))*sampDur;
            sampsToSecsMap(chunkIndex) = timeStamps;
            currentSample = chunkIndex(end)+1;
        end
          
        %% extract strobed events:
        % read the strobed word info (values & time stamps):
         switch rawFileType
            case 'plx'
                [~, strobedEvents.eventInfo.Ts, strobedEvents.eventInfo.Strobed] = plx_event_ts(rawFullPath, 257);
                % no start/stop...
            case 'pl2'
                strobedEvents.eventInfo = PL2EventTs(rawFullPath, 'Strobed');

                % read the time-stamps of recording start / stop events:
                strobedEvents.startTs = PL2StartStopTs(rawFullPath, 'start');
                strobedEvents.stopTs = PL2StartStopTs(rawFullPath, 'stop');

         end
        
        
    otherwise
        error('bad filetype. Time to reconsider your life choices');
end

%% subtract mean across channels:
if opts.commonAverageReferencing
    disp('Performing common average subtraction...')
    samplesMean = int16(mean(samples));
    % might wanna save out the samplesMean in case we want to view
    % it...
    % well, why not plot it:
    figure,
    plot(samplesMean);
    title('the samples mean, subtracted from all channels')
    % subtract:
    samples = bsxfun(@minus, samples, samplesMean);
end


%% remove artifacts:

if opts.removeArtifacts
    disp('Removing artifacts...')
    % set the standard deviation threshold:
    sdThresh            = 3.5;
    medAbs              = median(abs(single(samples)));
    sdMedAbs            = std(medAbs);
    if isfield(opts, 'removeArtifactsVisualize') && opts.removeArtifactsVisualize
        figure, 
        hold on
        plot(medAbs(1:1e2:end));
        hL(1) = line(xlim, [sdThresh*sdMedAbs sdThresh*sdMedAbs], 'Color', 'k');
        hL(2) = line(xlim, [sdThresh/2*sdMedAbs sdThresh/2*sdMedAbs], 'Color', 'r');
        hL(3) = line(xlim, [sdThresh*2*sdMedAbs sdThresh*2*sdMedAbs], 'Color', 'g');
        legend(hL, {'sdTh', 'sdTh/2', 'sdTh*2'})
        
    end
    % get artifact indices:
    idxBad              = median(abs(single(samples)) > (sdThresh*sdMedAbs));
    % remove the "bad" samples from 'samples' matrix:
    samples(:, idxBad)  = [];
    % remove the "bad" samples from timing vector too:
%     tsMap(idxBad)       = [];
    sampsToSecsMap(idxBad) = []; 
         
    fprintf('removed %0.0d of %0.0d samples, (%0.3f percent)\n', sum(idxBad), numel(idxBad), mean(idxBad)*1e2);
    
    if isfield(opts, 'plotProbeVoltage') && opts.plotProbeVoltage
        plot_probeVoltage(samples, 1e3);
        supertitle('AFTER ARTIFACT REMOVAL')
    end
        
end



%% Pack up and save:

% meta info:
info.dsn            = dsn;
info.rawFolder      = rawFolder;
info.rawFile        = rawFileName;
info.rawFullPath    = rawFullPath;
info.rawFileType    = rawFileType;
info.spkChNumber    = spkChNumber;
% info.strbChNumber   = strbChNumber;
info.opts           = opts;
info.datestr        = datestr(now, 'yyyymmddTHHMM');


% save info:
save(fullfile(opts.outputFolder, 'convertInfo.mat'),  'info');


% % timing data to mat file:
% disp('Saving mat file with timestamps & info')
% save(tsPath, 'sampsToSecsMap', 'info');

% save sampsToSecsMap (has to be 7.3 cause these can get BIG):
save(fullfile(opts.outputFolder, 'sampsToSecsMap.mat'),  'sampsToSecsMap', '-v7.3')

% save strobe info:
save(fullfile(opts.outputFolder, 'strobedEvents.mat'),  'strobedEvents')

% % ephys data to dat file:
fidout = fopen(datPath, 'a'); % opening file for appending
fwrite(fidout, samples, 'int16');
fclose(fidout);

fprintf('%f0.1s: CONVERSION COMPLETE!', toc)

dbclear if error

end
%% TEST ZONE
% clear t
% for iCh = 1:nChannels
%     tic;
%     pl = readPLXFileC(fullPathPlx, 'continuous', continuousChannelNumbers(iCh));
%     t.singleChannelLoad_readPlx(iCh) = toc;
%
%     tic;
%     [~,~,~,~, ad] = plx_ad(fullPathPlx, continuousChannelNumbers(iCh));
%     t.singleChannelLoad_plx_ad(iCh) = toc;
% end

%% look at activity
% nSecs = 3;
% fs = 40000;
% figure, hold on
% for iCh = 1:size(samples,2)
%     plot(iCh*500 + samples(1:10:nSecs*fs, iCh))
% end
% set(gca, 'XTick', 0:fs/2:nSecs, 'XTickLabel', 0:.5:nSecs)

% toc

%%

function out = adReadWrapper(fileName, chStr)

[~, ~, ~, ~, ad] = plx_ad(fileName, chStr);
out = int16(ad);

end


