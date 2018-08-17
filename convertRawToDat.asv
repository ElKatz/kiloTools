function [samples, datPath] = convertRawToDat(rawFullPath, opts)
%   [samples, datPath] = convertRawToDat(rawFullPath, opts)
%
% Converts continuous data from raw ephys file to matrix 'samples' of size:
% [nSamples, nChannels] and saves as binary .dat file in same folder as
% ephys file (unless specified otherwise in opts)
% INPUT:
%   ephysFullPath - Optional. full path to raw ephys data file.
%                   This file can be 'plx', 'mpx', 'oe', whatever, as long
%                   as it contains the continuous voltage traces.
%
% !!!!! AT THE MOMENT, THIS FUNCTION ONLY DEALS WITH plx FILES !!!!!
%r
%   opts - (optional) struct of options:
%   .outputFolder - full path to folder to save dat file (default: same
%                   folder as raw file)
%
% OUTPUT:
%   samples - [nChannels, nSamples] consisting of all continuous data
%   datPath - path to the .dat file


% 2do:
% generalize the identificaiton of continuous channel strings so that it
% works on plexon, alphaLab, openEphys etc...
%
% generalize to multiple file formats. vet.


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

% TIMESTAMPS & INFO: mat filename:
tsPath = fullfile(opts.outputFolder, 'timestampsAndInfo.mat');

%% begin conversion:

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
        %%% plx plx plx plx plx plx plx plx plx plx plx plx plx plx %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % create a list of all ad continuous channel names in cell array:
        [nCh, adChName]     = plx_adchan_names(rawFullPath);
        chNameList = cell(nCh,1);
        for ii = 1:nCh
            chNameList{ii} = adChName(ii,:);
        end
        
        % get indices for the spike channels (i.e. not the lfp):
        % *alphaLab spike files are saved with the name "SPK"  and plexon are
        % "CSPK" so I'm just gonna use "SPK" as the identifier:
        spkChannelStr = 'SPK';
        idxSpkCh      = false(numel(chNameList),1);
        for iCh = 1:numel(chNameList)
            if strfind(chNameList{iCh}, spkChannelStr)
                idxSpkCh(iCh) = true;
            else
                idxSpkCh(iCh) = false;
            end
        end
        % get number of spikes counts per ad channel and get indices for those that
        % have data:
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
        for iCh = 1:nChannels
            [~, n, ~, ~, ad] = plx_ad(rawFullPath, spkChNumber(iCh)); % returns signal in miliVolts
            if n>0
                tChRead(iCh) = toc;
                fprintf('\t%0.1fs: read channel #%0.0d \n', tChRead(iCh), spkChNumber(iCh));
                % data matrix 'samples':
                samples(iCh,1:n) = int16(ad);
            else
                warning('wtf?! something''s wrong i the assigning SPK channel numbers. fix this')
                keyboard
            end
        end
        
        %%  extract timing information from raw file in "real" time
        % we'll get:
        %   tsMap - time vector
        %   evTs - event timestamp (i.e. strobes)
        %   evSv - event strobe value
        %
        % 'tsMap' has a timestamp for every sample recorded. This will be a
        % vecotor of size nSamples. tsMap is used to convert from the spike 
        % index output of kiloSort (these are simply integers that indicate 
        % which sample number each spike occurred at), to time (in seconds)
        % relative to the beginning of the ephys recording.
        % This is needed because the event time stamps (evTs) from the raw
        % file are in same relative time (also in seconds).
        
        % get timestamps start values (tsStartVals) at start of each fragment:
        disp('Getting plexon timestamps for ad samples');
        [adfreq, ~, ts, fn] = plx_ad_gap_info(rawFullPath, spkChNumber(1));
        
        % loop over recording fragments
        tsMap = nan(sum(fn),1);
        currSample = 1;
        for ii = 1:length(fn)
            timeStamps = (0:(fn(ii)-1)) / adfreq;
            tsMap(currSample:(currSample - 1 + fn(ii))) = timeStamps + ts(ii);
            currSample = currSample + fn(ii);
        end
         
        disp('Getting plexon timestamps for strobed events');
        strbChNumber = 257; % !!! this is true for rig A/B. verify this is true for your rig too...
        [~, evTs, evSv] = plx_event_ts(rawFullPath, strbChNumber);
        
  
    otherwise
        error('bad filetype. Time to reconsider your life choices');
end

%% subtract mean across channels:
if opts.commonAverageReferencing
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
    % set the standard deviation threshold:
    sdThresh            = 10;
    sdMedAbs            = std(median(abs(single(samples))));
    idxBad              = median(abs(single(samples)) > (sdThresh*sdMedAbs));
    % remove the "bad" samples from 'samples' matrix:
    samples(:, idxBad)  = [];
    % remove the "bad" samples from timing vector too:
    tsMap(idxBad)       = [];
         
    fprintf('removed %0.0d of %0.0d samples, (%0.3f percent)\n', sum(idxBad), numel(idxBad), mean(idxBad)*1e2);
end



%% Pack up and save:

% meta info:
info.dsn            = dsn;
info.rawFolder      = rawFolder;
info.rawFile        = rawFileName;
info.rawFullPath    = rawFullPath;
info.rawFileType    = rawFileType;
info.spkChNumber    = spkChNumber;
info.strbChNumber   = strbChNumber;
info.opts           = opts;
info.datestr        = datestr(now, 'yyyymmddTHHMM');

% timing data to mat file:
disp('Saving mat file with timestamps & info')
save(tsPath, 'tsMap', 'evTs', 'evSv', 'info', '-v7.3');

% ephys data to dat file:
fidout = fopen(datPath, 'a'); % opening file for appending
fwrite(fidout, samples, 'int16');
fclose(fidout);

fprintf('%f0.1s: CONVERSION COMPLETE!', toc)

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