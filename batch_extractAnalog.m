%% batch_extractAnalog
%
% ad hoc script to go over datasets that have been converted and sorted
% before I had implemented the analog data extraction that is now part of
% the conversion process.

%% add necessary paths to toolboxes:
cd('D:\Code\Toolboxes\kiloTools');
paths = addPathsForSpikeSorting;
%% boolleans:

performConversion   = false;
performKiloSort     = true;

%% list of paths to raw ephys files

fs              = 40000;
nCh             = 32;
Nfilt           = 32*3;     % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)

% set your options:
probeGeometry                   = 'linear50'; % see probeGeometry2coords.m for options
opts.commonAverageReferencing   = false;
opts.specificChannels           = false; % either false or (eg) 65:96
if opts.specificChannels
    nCh = length(opts.specificChannels);
end
opts.plotProbeVoltage           = true;
opts.removeArtifacts            = false;

% input your folders:
folderList = {...
%   'Y:\LAB PROJECTS\scInactivateAndRecord\data\ramsey20190920\'...
    'Y:\LAB PROJECTS\scInactivateAndRecord\data\ramsey20200123\'...
%     'Y:\LAB PROJECTS\scInactivateAndRecord\data\ramsey20200121\'...
%     'Y:\LAB PROJECTS\scInactivateAndRecord\data\ramsey20200122\'...
    };
nFiles = numel(folderList);

% Make list of dat files by adding .dat and inserting 'kiloSorted' folder:
rawPath         = cell(nFiles,1);
datPathList     = cell(nFiles,1);
kiloFolderList  = cell(nFiles,1);
fileName        = cell(nFiles,1);
for iF = 1:nFiles
    fileList    = dir(folderList{iF});
    idxGood     = arrayfun(@(x) ~strcmp(x.name(1), '.'), fileList); % identify files that do not start with '.'
    fileList    = fileList(idxGood);
    idxPl2      = arrayfun(@(x) ~isempty(strfind(x.name, '.pl2')), fileList); % identify file that contain .pl (i.e. .plx or .pl2) 
    rawFile     = fileList(idxPl2).name;
    rawPath{iF} = fullfile(folderList{iF}, rawFile);
    [~, fileName{iF}]   = fileparts(rawPath{iF});
    kiloFolderList{iF}  = fullfile(folderList{iF}, 'kiloSorted');
    datPathList{iF}     = fullfile(kiloFolderList{iF}, [fileName{iF} '.dat']);
end

%%  (1) convert, (2) copy masterMegaFile into kiloSorted folder, (3) sort:

for iF = 1:nFiles
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(['~~~~~  ' rawPath{iF}])
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    % convert:
    if performConversion
        if exist(rawPath{iF}, 'file')
            convertRawToDat(rawPath{iF}, opts);
        end
    end
    
    % kiloSort:
    if performKiloSort
        copyfile(fullfile(paths.kiloTools, 'masterMegaFile.m') ,kiloFolderList{iF});
        if exist(datPathList{iF}, 'file')
            cd(kiloFolderList{iF})
            masterMegaFile(datPathList{iF}, fs, nCh, probeGeometry, Nfilt);
        end
    end
    
end

try
    load gong
    soundsc(y, Fs)
%     evalc('soundsc(double(aiffread(''/System/Library/Sounds/Glass.aiff''))'',50000)');
catch
end

%% NOW YOU SORT BY HAND.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% manually curate spikes with PHY %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% and then run batch_mkKsFigs.m %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% fin.


%%%%%%%%%%%%%%%%%%%%%


