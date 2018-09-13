function [] = extractTimingAndStrobes(rawFullPath)

% HACKKKKK!!!!!


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

%%


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
        
        % must read in a spike channel to construct the "timestamp map" from
        % samples (kilosort) to time in seconds.
        ad = PL2Ad(rawFullPath, 'SPKC01');

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
          
        % read the strobed word info (values & time stamps).
        strobedEvents.eventInfo = PL2EventTs(rawFullPath, 'Strobed');
        strobedEvents.RSTARTInfo = PL2EventTs(rawFullPath, 'RSTART');
        strobedEvents.RSTOPInfo = PL2EventTs(rawFullPath, 'RSTOP');
        
        % save strobedEvents:
save(fullfile(opts.outputFolder, 'strobedEvents.mat'),  'strobedEvents');
% save sampsToSecsMap (has to be 7.3 cause these can get BIG):
save(fullfile(opts.outputFolder, 'sampsToSecsMap.mat'),  'sampsToSecsMap', '-v7.3')