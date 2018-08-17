function [datPath] = convertToDat_plxFolder(plxFolderPath, opts)
%   [datPath] = convertToDat_plxFolder(plxFullPath, opts)
%
% This function is similar to its predecessor 'convertToDat_plxFile' but 
% instead of operating using a single plx file that has all data as input, 
% it takes a folder that has seperate files for each channels. 
% *** this is critical. the function only works on file that have been
% split by plxUtil into 1 file per channel ***
%
%   SPKC    - continuous spike channel (ie highpassed wideband)
%   FP      - field potential (ie lowpassed wideband) 
%   AN      - analog input channels (eg eyeX, joystick)
%   RSTART  - all start times (that initiate ephys recording
%   RSTOP   - all stop times (that pause ephys recording)
%   STROBED - all strobe times and values.
% it loads SPKC and converts into a [nSamples, nChannels] matrix and saves 
% as binary .dat file in same folder under "kiloSroted"
% INPUT:
%   plxFolderPath - Optional. full path to raw ephys plx file.
%                   Must contain the continuous voltage traces.
%
%   opts - (optional) struct of options:
%
% OUTPUT:
%   datPath - path to the output .dat file


%% paths:
addPathsForSpikeSorting;

%% data file & folder names:

if ~exist('plxFullPath', 'var')
    [plxFolder] = uigetdir('Select folder for conversion', '~/Dropbox/Code/spike_sorting/');
end

[~, dsn] = fileparts(plxFolder);

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
% all continuous voltage data (SPKC) gets saved in the binary .dat file.
% all strobed event timestamps ('ts') get saved in a .mat file.
% meta info also gets saved there too. both are named after the 'dsn':

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

% get file list and index into those that have the dsn (taking only difrst
% 11 chars which should be of format subYYYYMMDD):
fileList        = dir(plxFolder);
idx             = arrayfun(@(x) ~isempty(strfind(x.name, dsn(1:11))), fileList);
fileName    = {fileList(idx).name};

plxFiles = struct(...
    'fileName', [], ...
    'chName', [], ...
    'isSpkc', [], ...
    'chNumber', [], ...
    'nSamples', []);


for iF = 1:numel(fileName)
    
    % full path to plx file:
    plxFullPath = fullfile(plxFolder, fileName{iF});
    
    % file name:
    plxFiles(iF).fileName = fileName{iF};

    if strfind(fileName{iF}, 'SPKC')
        plxFiles(iF).isSpkc = true;
        
        % find the channel that has samples (should be a single channel):
        [~, samplecounts] = plx_adchan_samplecounts(plxFullPath);
        chPtr = find(samplecounts~=0);
        assert(numel(chPtr)==1, 'oh oh- more than one channel had data. that bad')
        plxFiles(iF).nSamples = samplecounts(chPtr);
        
        % get channel name: 
        [~, names] = plx_adchan_names(plxFullPath);
        plxFiles(iF).chName = names(chPtr,:);
        
        % get channel plexon-based number:
        [~, numbers] = plx_ad_chanmap(plxFullPath);
        plxFiles(iF).chNumber = numbers(chPtr);
    else
        plxFiles(iF).isSpkc = false;
        
        if strfind(fileName{iF}, 'RSTART')
            plxFiles(iF).chName = 'Start';
            plxFiles(iF).chNumber = plx_event_resolve_channel(plxFullPath, 'Start');
            
        elseif strfind(fileName{iF}, 'RSTOP')
            plxFiles(iF).chName = 'Stop';
            plxFiles(iF).chNumber = plx_event_resolve_channel(plxFullPath, 'Stop');
            
        elseif strfind(fileName{iF}, 'Strobed')
            plxFiles(iF).chName = 'Strobed';
            plxFiles(iF).chNumber = plx_event_resolve_channel(plxFullPath, 'Strobed');
        end
    end
end

%% build data matrix 'samples' of size [nChannels, nSamples]:

% get all SPKC channels:
idxSpkc     = [plxFiles.isSpkc];
nChannels   = sum(idxSpkc);
assert(numel(unique([plxFiles(idxSpkc).nSamples]))==1, 'oh oh, different SPKC channels have a different number of samples!');
nSamples = unique([plxFiles(idxSpkc).nSamples]);

% prealloc:
samples     = zeros(nChannels, nSamples, 'int16');
tChRead     = nan(nChannels,1); % time keeping

disp(['Getting data from ' num2str(nChannels) ' spike channels!'])
spkcFileNames = {plxFiles(idxSpkc).fileName};
spkcChNumbers = [plxFiles(idxSpkc).chNumber];
for iCh = 1:nChannels
    tic
    plxFullPath = fullfile(plxFolder, spkcFileNames{iCh});
    [~, n, ~, ~, ad] = plx_ad(plxFullPath, spkcChNumbers(iCh)); % returns signal in miliVolts
    if n>0
        tChRead(iCh) = toc;
        fprintf('\t%0.3f secs to read channel #%0.0d \n', tChRead(iCh), spkcChNumbers(iCh));
        % data matrix 'samples':
        samples(iCh,1:n) = int16(ad);
    else
        warning('wtf?! something''s wrong in the assigning SPK channel numbers. fix this')
        keyboard
    end
end

%% remove artifacts:

if opts.removeArtifacts
    sdThresh = 12;
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
[adfreq, ~, ts, fn] = plx_ad_gap_info(plxFullPath, spkcChNumbers(1));

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
idxStrobed = arrayfun(@(x) ~isempty(strfind(x.fileName, 'Strobed')), plxFiles);
[~, evTs, evSv] = plx_event_ts(plxFullPath, plxFiles(idxStrobed).chNumber);


%% Pack up and save:
% add meta info:
info.dsn            = dsn;
info.rawFolder      = plxFolder;
info.plxFiles       = plxfiles;
info.opts           = opts;
info.datestr        = datestr(now, 'yyyymmddTHHMM');

% save timing data to mat file:
disp('Saving mat file with timestamps & info')
save(tsPath, 'tsMap', 'evTs', 'evSv', 'info', '-v7.3');

% save ephys data to .dat file:
fidout = fopen(datPath, 'w'); % opening file for appending
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