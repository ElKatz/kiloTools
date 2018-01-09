%% Define paths to your dataset:
dirDataset  = '~/Dropbox/_transfer_big_data/20170615/data/convertedToRaw/';
dsn         = 'f170615-0002'; % dataset name

% add paths:
addpath(genpath('~/Dropbox/Code/spike_sorting/packages/KiloSort')) % path to kilosort folder
addpath(genpath('~/Dropbox/Code/spike_sorting/toolboxes/npy-matlab')) % path to npy-matlab scripts

% run 'createChanMap.m':
pathToYourCreateChanMap = fullfile(dirDataset, 'createChannelMap.m'); % take from Github folder and put it somewhere else (together with the master_file)
run(pathToYourCreateChanMap)

% run 'standardConfig.m':
pathToYourConfigFile = fullfile(dirDataset, 'standardConfig.m'); % take from Github folder and put it somewhere else (together with the master_file)
run(pathToYourConfigFile)

%% sort that sweet sweet data:
if ops.GPU     
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end
tic;
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

%% save preAutoMerge sort results?
if isfield(ops, 'savePreAutoMerge') && ops.savePreAutoMerge
    % save the preAutoMerge results:
    preAutoMergeFolder = fullfile(dirDataset, 'preAutoMerge');
    mkdir(preAutoMergeFolder)
    save(fullfile(preAutoMergeFolder,  'rez.mat'), 'rez', '-v7.3');
    rezToPhy(rez, preAutoMergeFolder);
end

%% AutoMerge:
% rez2Phy will use clusters from the new 5th column of st3 if you run this
rez = merge_posthoc2(rez);

% save post autoMerge results in root folder:
save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3');
rezToPhy(rez, ops.root);

% remove temporary file
delete(ops.fproc);
ding
