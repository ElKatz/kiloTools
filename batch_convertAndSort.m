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
cd('D:\Code\Toolboxes\kiloTools');
paths = addPathsForSpikeSorting;
%% boolleans:

performConversion   = true;
performKiloSort     = true;

%% list of paths to raw ephys files

fs              = 40000;
nCh             = 28;
probeGeometry   = 'linear125';


folderList = {...
    'Z:\LAB PROJECTS\cla_record\data\\p20190323a\', ...
    'Z:\LAB PROJECTS\cla_record\data\\p20190323b\', ...
    };

nFiles = numel(folderList);



% Make list of dat files by adding .dat and inserting 'kiloSorted' folder:
rawPath         = cell(nFiles,1);
datPathList     = cell(nFiles,1);
kiloFolderList  = cell(nFiles,1);
fileName        = cell(nFiles,1);
for iF = 1:nFiles
    fileList    = dir(folderList{iF});
    idxPl2      = arrayfun(@(x) ~isempty(strfind(x.name, 'pl2')), fileList);
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
            convertRawToDat(rawPath{iF});
        end
    end
    
    % kiloSort:
    if performKiloSort
        copyfile(fullfile(paths.kiloTools, 'masterMegaFile.m') ,kiloFolderList{iF});
        if exist(datPathList{iF}, 'file')
            cd(kiloFolderList{iF})
            masterMegaFile(datPathList{iF}, fs, nCh, probeGeometry);
        end
    end
    
end

%% NOW YOU SORT BY HAND.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% manually curate spikes with PHY %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% and then run batch_mkKsFigs.m %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%
% % %% AND THEN:
% % 
% % %%
% % 
% % 
% % for iF = 1:nFiles
% % 
% %     % get strobes:
% %     strobedEvents = getStrobedFromRaw(rawPath{iF});
% %     % save strobedEvents:
% %     save(fullfile(kiloFolderList{iF}, 'strobedEvents.mat'),  'strobedEvents');
% % 
% %     
% %     sp = getSp(kiloFolderList{iF});
% %     
% %     % combo:
% %     combo.sp = sp;
% %     combo.strobedEvents = strobedEvents;
% %     save(fullfile(folderList{iF}, [fileName{iF} '_combo.mat']), '-struct', 'combo',  '-v7.3');
% % end
% %     
    %%

% Fuckup:
% I acccidentally resorted these files after having already sorted+phy'd
% them. The resorting (with no phy) overwrote most files. most notable,
% spike_clusters.npy (which Amar uses to read in phy'd clusters). So it is
% crucial that we read the .csv and and not the .npy for these sessions!
% {'Z:\katz_server\fstAttention\data\fst_with_sc_inactivation\20171214\pre\20171214_t1240' }
%     {'Z:\katz_server\fstAttention\data\fst_with_sc_inactivation\20171211\pre\20171211_t1253' }
%     {'Z:\katz_server\fstAttention\data\fst_with_sc_inactivation\20171211\post\20171211_t1655'}
