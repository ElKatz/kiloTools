function [raw, rawInfo] = convertToRaw_channelByChannel(dataFileType, opts)
%
% NEED TO UPDATE!!!!
%
% [raw, rawInfo] = convertToRaw(dataFileType, opts)
%
% Converts ephys files to "raw" matrix of [nChannels, nSamples] and saves
% the file as a binary .dat file (the preferred format for Kilosort).
%
% Function will prompt user to select files for conversion (either multiple
% files or one large one. *if you run into memory issues with large files
% consdier spilitting).
%
% Function default is to segment the files
% 3. identify start/stop TTL pulses and segment
% 4. append to dat



% fullPathPlx = '/Users/katzln/Dropbox/Code/spike_sorting/pilot_datasets/tst_mpx2plx_vs_mpx2mat/F170418-0002 - CopyPost.plx';
% BIGGEST 2DO:
% KEEP TRACK OF TIMING!!!!!
% SEGMENT TTL PULSE
% gotta delete any existing .dat files or increase iterator..



% mpx = load('~/Dropbox/Code/spike_sorting/pilot_datasets/leor_vProbe_attn_20170511_cont/f170511-0019.mat');

%% input:
if ~exist('dataFileType', 'var')
    fileTypeOptions = {'plx', 'matFromMpx'};
    disp('You did not supply a dataFileType sillyhead!')
    disp('Please select one of the following:')
    disp(fileTypeOptions)
    dataFileType = input('Enter your dataFileType now:\n');
end
% default opts strcut:
if ~exist('opts', 'var'),                   opts = struct;                      end
if ~isfield(opts, 'saveDat'),               opts.saveDat = true;                end
if ~isfield(opts, 'saveMat'),               opts.saveMat = false;                end
if ~isfield(opts, 'performSegmentation'),   opts.performSegmentation = false;    end

%% paths:
addpath(genpath('~/Dropbox/Code/spike_sorting/toolboxes/Matlab Offline Files SDK'));
addpath(genpath('~/Dropbox/Code/repos/pdstools/dependencies'));
addpath(genpath('~/Dropbox/Code/spike_sorting/toolboxes'));

%% data folder:

if ~isfield(opts, 'startFolder')
    opts.startFolder = '~/Dropbox/Code/spike_sorting/pilot_datasets/';
end

% get the data files:
switch dataFileType
    case 'plx'
        [dataFiles, dataFolder] = uigetfile('*.plx', 'Select file/s for conversion', opts.startFolder, 'MultiSelect', 'on');
        chString_spike = 'SPKC';
        chString_lfp   = 'FP';
        
    case 'matFromMpx'
        [dataFiles, dataFolder] = uigetfile('*.mat', 'Select file/s for conversion', opts.startFolder, 'MultiSelect', 'on');
        chString_spike = 'CSPK';
        chString_lfp   = 'CLFP';
        
    otherwise
        warning('invalid dataFileType arguemnt')
        disp('BOOOO!')
        disp('User must supply fileType. options are:')
        disp(fileTypeOptions)
        return;
end

% make sure dataFile is a cell:
if ~iscell(dataFiles)
    dataFiles = mat2cell(dataFiles, 1);
end
% sort by name:
dataFiles = sort(lower(dataFiles));

nFiles = numel(dataFiles);

% display:
disp('--------------------------------------------------------------')
disp('converting file/s from folder:')
disp(dataFolder)
disp('--------------------------------------------------------------')
disp(['Found ' num2str(nFiles) ' files:'])
disp(dataFiles')

dsn = dataFiles{1}(1:end-4); % datasetname

% raw binary .dat file name is named after the first plx file:
% spkDataFileName = [dsn '.dat'];
spkTimeFileName = [dsn '_times.dat'];


disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp(['FROM THE ASHES OF ' dataFileType '... LIKE THE PHOENIX... I GIVE YOU:'])
disp(['   '  dsn])
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

%% begin conversion:

tic;

% prealloc spkTime, event time stamps (ts) and strobe value (sv):
spk     = [];
tSpk    = [];
lfp     = [];
tLfp    = [];
tEvent  = [];
vEvent  = [];

nFiles = numel(dataFiles);
for iF = 1:nFiles
    % full path to file:
    fullPathData = fullfile(dataFolder, dataFiles{iF});
    fprintf('\n\n%s:  performing conversion', dataFiles{iF})
    
    %% %%%%% GET THAT SWEET SWEET DATA %%%%%%%%%
    % %
    % %     % extract raw data from spike channels:
    % %     [rawSpk] = extractRawDataFromFile(fullPathData, dataFileType, chString_spike, opts);
    % %
    % %     % extract raw data from LFP channels:
    % %     [rawLfp] = extractRawDataFromFile(fullPathData, dataFileType, chString_lfp, opts);
    
    mpx     = load(fullPathData);
    nChannels = 24;
    
      % enforce equal sampling rate accross channels:
        mpx = enforceEqualSmplingRateAcrossChannels(mpx, 'SPKC');
        
    %% code depends on whether mpx-to-mat conversion was in struct or not:
    if isstruct(mpx.CSPK_001)
        nSamples = numel(mpx.CSPK_001.Samples);
        samples = zeros(nChannels, nSamples, 'int16');
        for iCh = 1:nChannels
            samples(iCh,:) = mpx.(['CSPK_0' sprintf('%.2d', iCh)]).Samples;
        end
        tSamples = linspace(mpx.CSPK_001.TimeBegin, mpx.CSPK_001.TimeEnd, nSamples);
        % get event timestmpas (ts) and strobe values (sv):
        if isfield(mpx, 'CInPort_001')
            tEventFile  = mpx.CInPort_001.Samples(1,:) ./ (mpx.CInPort_001.KHz*1e3);
            vEventFile  = mpx.CInPort_001.Samples(2,:);
        else
            tEventFile = [];
            vEventFile = [];
        end
    else
        nSamples = numel(mpx.CSPK_001);
        samples = zeros(nChannels, nSamples, 'int16');
        for iCh = 1:nChannels
            samples(iCh,:) = mpx.(['CSPK_0' sprintf('%.2d', iCh)]);
        end
        tSamples = linspace(mpx.CSPK_001_TimeBegin, mpx.CSPK_001_TimeEnd, nSamples);
        % get event timestmpas (ts) and strobe values (sv):
        if isfield(mpx, 'CInPort_001')
            tEventFile  = mpx.CInPort_001(1,:) ./ (mpx.CInPort_001_KHz*1e3);
            vEventFile  = mpx.CInPort_001(2,:);
        else
            tEventFile = [];
            vEventFile = [];
        end
    end
    
    
    %% appending data to files on harddrive (only performed from spkData and spkTime:)
    if isfield(opts, 'saveDat') && opts.saveDat
        fprintf('\nAppending to file...')
        
        
        for iCh = 1:nChannels
            % spkData:
            spkDataFileName = [dsn '-ch' sprintf('%.2d', iCh) '.dat'];
            spkDataPath = fullfile(dataFolder, spkDataFileName);
            fidout = fopen(spkDataPath, 'a'); % opening file for appending
            fwrite(fidout, samples(iCh,:), 'int16');
            fclose(fidout);
        end
        
        % spkTime:
        spkTimePath = fullfile(dataFolder, spkTimeFileName);
        fidout = fopen(spkTimePath, 'a'); % opening file for appending
        fwrite(fidout, tSamples, 'double');
        fclose(fidout);
        fprintf('\tdone');
        
    end
    
    %% concatenting data from each file:
    
    if nargout > 0 || (isfield(opts, 'saveMat') && opts.saveMat)
        
        % values:
        %         spk = [spk, rawSpk.samples];
        %         lfp = [lfp, rawLfp.samples];
        
        %         % times:
        %         tSpk = [tSpk(:); rawSpk.tSamples(:)];
        %         tLfp = [tLfp(:); rawLfp.tSamples(:)];
        %
        % event time & value:
        tEvent = [tEvent(:); tEventFile(:)];
        vEvent = [vEvent(:); vEventFile(:)];
    end
    
    fprintf('\nTime = %0.1f, file %d of %d complete', toc, iF, nFiles)
    a =whos; 
    totalGB = sum(arrayfun(@(x) x.bytes, a)) / 1024 /1024 / 1024;
    fprintf('\nTotal GB: %d', totalGB);
end % iF
fprintf('\nTime = %0.1f, Conversion complete!', toc)

%% save/output

if nargout > 0 || (isfield(opts, 'saveMat') && opts.saveMat)
    
    % get all meta info:
    rawInfo.dsn             = dsn;
    rawInfo.dataFolder      = dataFolder;
    rawInfo.dataFiles       = dataFiles;
    rawInfo.dataType        = dataFileType;
    %     rawInfo.chNameSpk       = rawSpk.chNameList;
    %     rawInfo.chNameLfp       = rawLfp.chNameList;
    rawInfo.spkDataPath     = spkDataPath;
    rawInfo.opts            = opts;
    rawInfo.datestr         = datestr(now, 'yyyymmddTHHMM');
    
    % get all data:
    raw.spk     = spk;
    raw.tSpk    = tSpk;
    raw.lfp     = lfp;
    raw.tLfp    = tLfp;
    raw.tEvent  = tEvent;
    raw.vEvent  = vEvent;
    raw.info    = rawInfo;
    
    % save:
    if (isfield(opts, 'saveMat') && opts.saveMat)
        fprintf('\nsaving mat files...')
        save(fullfile(dataFolder, 'raw'), 'raw', '-v7.3');
        fprintf('\nTime = %0.1f, Saving complete!', toc)
    end
end

fprintf('\nDONE AND DONE!\n\n')
ding
%% vis
% figure,
% plot(raw.tSpk(1:1e4:end), raw.spk(1,1:1e4:end))




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
% %
% %     % identify the continuous channels by name:
% %     [nCh, adChName]     = plx_adchan_names(fullPathPlx);
% %     [~,   adChNumber]   = plx_ad_chanmap(fullPathPlx);
% %     chNameList = cell(nCh,1);
% %     for ii = 1:nCh
% %         chNameList{ii} = adChName(ii,:);
% %     end
% %     % continuous spike channels are identified by the strong below. This is
% %     % not a great method but it works for alphaLab files (which are named
% %     % "SPK") and plexon recorded files (named "CSPK")
% %     continuousChannelStr = 'SPK';
% %     % get indices for continuous channels:
% %     idx         = false(numel(chNameList),1);
% %     for iCh = 1:numel(chNameList)
% %         if strfind(chNameList{iCh}, continuousChannelStr)
% %             idx(iCh) = true;
% %         else
% %             idx(iCh) = false;
% %         end
% %     end
% %     contChNumber = adChNumber(idx); % adding 1 because adChNumber is zero-based list
% %
% %     % make sure that the channel numbers are similar between files:
% %     if iF==1
% %         contChNumberFirstFile = contChNumber;
% %     else
% %         assert(mean(contChNumber==contChNumberFirstFile)==1, 'goddamn! There''s a mismatch between the continuous channel numbers of this file vs 1st file. wtf?');
% %     end
% %
% %     % build data matrix 'samples' of size [nChannels, nSamples]:
% %     nChannels   = numel(contChNumber);
% %     [~, ~, ~, contCounts] = plx_info(fullPathPlx, 1);
% %     nSamples    = contCounts(contChNumber(1));
% %     samples     = zeros(nChannels, nSamples, 'int16');
% %     tChRead     = nan(nChannels,1);
% %     chGood      = false(nChannels,1);
% %     hWait = waitbar(0, 'Loading up channel by channel and populating ''samples''');
% %     for iCh = 1:nChannels
% %         tic
% %         [adfreq, n, ts, fn, ad] = plx_ad(fullPathPlx, contChNumber(iCh)); % retunrs signal in miliVolts
% %         if n>0
% %             chGood(iCh) = true;
% %             tChRead(iCh) = toc;
% %             sprintf('Took me %0.3f secs to read channel %0.0d', tChRead(iCh), contChNumber(iCh));
% %
% %             % data matrix 'samples':
% %             samples(iCh,1:n) = int16(ad);
% %         else
% %             chGood(iCh) = false;
% %         end
% %         waitbar(iCh/nChannels, hWait)
% %     end
% %     close(hWait)
% %     disp(['Got data from' num2str(sum(chGood)) ' continuous channels!'])
% %
% %     % take only good channels:
% %     samples = samples(chGood, :);




