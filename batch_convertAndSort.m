%% batch_convertAndSort
%
% example script for batch conversion and sorting of raw ephys files. 
%
% INSTRUCTIONS: 
% - provide script with list of files.
% script will convert each raw file into .dat, will copy masterMegaFile.m 
% into each directory, and kiloSort each
% - make sure you have correct path to kiloTools

%% add necessary paths to toolboxes:
paths = addPathsForSpikeSorting;

%% list of paths to files 
% do not add file extension, e.g. 'C:\myFolder\myFile'

pathList = {...
    'C:\EPHYS\Data\Katz\test\20171211_t2055666', ...
    'C:\EPHYS\Data\Katz\test\20171211_t2055', ...
    
    };

nFiles = numel(pathList);

%% extensions are added here:
rawPathList = cellfun(@(x) [x '.plx'], pathList, 'UniformOutput', 0);
datPathList = cellfun(@(x) [x '.dat'], pathList, 'UniformOutput', 0);

%% copy masterM<egaFile into each folder:
folderList = cellfun(@(x) fileparts(x), pathList, 'UniformOutput', 0);
for iF = 1:nFiles
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(['~~~~~~~  ' pathList{iF} '  ~~~~~~~'])
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    if exist(folderList{iF}, 'dir')
        % convert:
        if exist(rawPathList{iF}, 'file')
            convertRawToDat(rawPathList{iF});
        else
            warning(['Could not find raw ephys file: ' rawPathList{iF}])
            continue;
        end
        
        % copy masterMegaFile:
        copyfile(fullfile(paths.kiloTools, 'masterMegaFile.m') ,folderList{iF});
        
        % kiloSort:
        if exist(datPathList{iF}, 'file')
            masterMegaFile(datPathList{iF});
        else
            warning(['Could not find dat  file: ' datPathList{iF}])
            continue;
        end
    else
        warning(['Could not find folder: ' folderList{iF}])
    end
end
%% convert & sort

