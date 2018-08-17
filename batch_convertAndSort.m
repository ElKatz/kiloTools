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
%% boolleans:

performConversoin   = true;
performKiloSort     = true;
%% list of paths to raw ephys files

% 'Z:\LAB PROJECTS\fstAttention\data\fst\s20170711a-non-non\20170711_1_2700sec_onwards_where_stimuli_are_in_correct_location.plx',...
% 'Z:\LAB PROJECTS\fstAttention\data\fst\s20171030a-non-non\20171030_t1351.plx',...
rawPathList = {...
    
%     'D:\Data\katz\GPe_recording\20180808_pilot_vProbe_with_PLDAPS_vK2\20180808_pilot_vProbe_with_PLDAPS_vK2.pl2',...
%     'D:\Data\katz\GPe_recording\20180808_pilot_vProbe_with_PLDAPS_vK2_withStartAndStop\20180808_pilot_vProbe_with_PLDAPS_vK2_withStartAndStop.pl2',...
        'D:\Data\katz\GPe_recording\20180814a\20180814_t1407.pl2', ...
        'D:\Data\katz\GPe_recording\20180814b\20180814_t1421.pl2', ...
    }; 

nFiles = numel(rawPathList);

fs              = 40000;
nCh             = 31;
probeGeometry   = 'linear200';

%% Make list of dat files by adding .dat and inserting 'kiloSorted' folder:
datPathList     = cell(nFiles,1);
kiloFolderList  = cell(nFiles,1);
for iF = 1:nFiles
    [folder, file] = fileparts(rawPathList{iF});
    kiloFolderList{iF}  = fullfile(folder, 'kiloSorted');
    datPathList{iF}     = fullfile(kiloFolderList{iF}, [file '.dat']);
end

%% (1) convert, (2) copy masterMegaFile into kiloSorted folder, (3) sort:

for iF = 1:nFiles
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(['~~~~~  ' rawPathList{iF}])
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    % convert:
    if performConversoin
        if exist(rawPathList{iF}, 'file')
%             convertRawToDat(rawPathList{iF});
            [samples] = convertRawToDat_beforeChange(rawPathList{iF});
        else
            warning(['Conversion fail. Could not find raw ephys file: ' rawPathList{iF}])
            continue;
        end
    end
    
    % kiloSort:
    if performKiloSort
        copyfile(fullfile(paths.kiloTools, 'masterMegaFile.m') ,kiloFolderList{iF});
        if exist(datPathList{iF}, 'file')
            cd(kiloFolderList{iF})
            masterMegaFile(datPathList{iF}, fs, nCh, probeGeometry);
        else
            warning(['kiloSort fail. Could not find dat file: ' datPathList{iF}])
            continue;
        end
    end
    
end


%%

% Fuckup:
% I acccidentally resorted these files after having already sorted+phy'd
% them. The resorting (with no phy) overwrote most files. most notable,
% spike_clusters.npy (which Amar uses to read in phy'd clusters). So it is
% crucial that we read the .csv and and not the .npy for these sessions!
% {'Z:\katz_server\fstAttention\data\fst_with_sc_inactivation\20171214\pre\20171214_t1240' }
%     {'Z:\katz_server\fstAttention\data\fst_with_sc_inactivation\20171211\pre\20171211_t1253' }
%     {'Z:\katz_server\fstAttention\data\fst_with_sc_inactivation\20171211\post\20171211_t1655'}
