function [datPath] = convertPlxToDat(plxFullPath, opts)
%   [datPath] = convertPlxToDat(plxFullPath, opts)
%
% Converts continuous data from plx ephys file to matrix 'samples' of size:
% [nSamples, nChannels] and saves as binary .dat file in same folder as
% ephys file (unless specified otherwise in opts)
% INPUT:
%   plxFullPath - Optional. full path to raw ephys plx file.
%                 Must contain the continuous voltage traces.
%
%   opts - (optional) struct of options:
%
% OUTPUT:
%   datPath - path to the .dat file


%% paths:
addPathsForSpikeSorting;

%% data file & folder names:

if ~exist('plxFullPath', 'var')
    [plxFileName, plxFolder] = uigetfile('*.*', 'Select plx/pl2 files for conversion', '~/Dropbox/Code/spike_sorting/');
else
    [plxFolder, plxFileName, plxFileType]  = fileparts(plxFullPath);
    plxFileName = [plxFileName plxFileType]; % have name include extension
    plxFileType = plxFileType(2:end); % remove dot
end

% full path to plx file:
plxFullPath = fullfile(plxFolder, plxFileName);

% datasetname:
dsn = plxFileName(1:end-4); % name of plx file without extension

%% use optional arguements and/or set defaults:
% init:
if ~exist('opts', 'var')
    opts = struct;
end

% output folder called 'kiloSorted' goes inside current folder:
outputFolder = fullfile(plxFolder, 'kiloSorted');
mkdir(outputFolder);

% remove artifacts
if ~isfield(opts, 'removeArtifacts')
    opts.removeArtifacts = true;
end

%% .dat & .mat:
% all ephys data gets saved in the binary .dat file.
% all strobed event timestamps ('ts') get saved in a .mat file.
% meta info also gets saved there. both are named after the 'dsn':

% .dat:
datPath = fullfile(outputFolder, [dsn '.dat']);
if exist(datPath, 'file')
    delete(datPath);
end

% .mat:
tsPath = fullfile(outputFolder, [dsn '.mat']);

%% begin conversion:

disp('--------------------------------------------------------------')
fprintf('Performing conversion of %s\n', dsn)
disp('--------------------------------------------------------------')

% create a list of all ad continuous channel names in cell array:
[nCh, adChName]     = plx_adchan_names(plxFullPath);
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
[~, samplecounts] = plx_adchan_samplecounts(plxFullPath);
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
[~,   adChNumber]   = plx_ad_chanmap(plxFullPath);
spkChNumber = adChNumber(idxGoodCh);

disp(['Getting data from ' num2str(sum(idxGoodCh)) ' spike channels!'])
for iCh = 1:nChannels
    tic
    [~, n, ~, ~, ad] = plx_ad(plxFullPath, spkChNumber(iCh)); % returns signal in miliVolts
    if n>0
        tChRead(iCh) = toc;
        fprintf('\tTook me %0.3f secs to read channel #%0.0d \n', tChRead(iCh), spkChNumber(iCh));
        % data matrix 'samples':
        samples(iCh,1:n) = int16(ad);
    else
        warning('wtf?! something''s wrong i the assigning SPK channel numbers. fix this')
        keyboard
    end
end

%% remove artifacts:

if opts.removeArtifacts
    sdThresh = 10;
    sdMedAbs = std(median(abs(single(samples))));
    idxBad   = median(abs(single(samples)) > (sdThresh*sdMedAbs));
    samples(:, idxBad) = [];
    fprintf('removed %0.0d of %0.0d samples, (%0.3f percent)\n', sum(idxBad), numel(idxBad), mean(idxBad)*1e2);
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
[adfreq, ~, ts, fn] = plx_ad_gap_info(plxFullPath, spkChNumber(1));

% loop over recording fragments
tsMap = nan(sum(fn),1);
currSample = 1;
for ii = 1:length(fn)
    timeStamps = (0:(fn(ii)-1)) / adfreq;
    tsMap(currSample:(currSample - 1 + fn(ii))) = timeStamps + ts(ii);
    currSample = currSample + fn(ii);
end

% remove artifacts samples from timing vector too:
if opts.removeArtifacts
    tsMap(idxBad) = [];
end

disp('Getting plexon timestamps for strobed events');
strbChNumber = 257; % !!! this is true for rig A/B. verify this is true for your rig too...
[~, evTs, evSv] = plx_event_ts(plxFullPath, strbChNumber);


%% Pack up and save:
% add meta info:
info.dsn            = dsn;
info.rawFolder      = plxFolder;
info.rawFile        = plxFileName;
info.rawFullPath    = plxFullPath;
info.rawFileType    = plxFileType;
info.spkChNumber    = spkChNumber;
info.strbChNumber   = strbChNumber;
info.opts           = opts;
info.datestr        = datestr(now, 'yyyymmddTHHMM');

% save timing data to mat file:
disp('Saving mat file with timestamps & info')
save(tsPath, 'tsMap', 'evTs', 'evSv', 'info', '-v7.3');

% save ephys data to .dat file:
fidout = fopen(datPath, 'a'); % opening file for appending
fwrite(fidout, samples, 'int16');
fclose(fidout);

disp('CONVERSION COMPLETE!')

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