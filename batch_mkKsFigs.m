%% batch_mkKsFigs
%
% make a bunch of figures for each kiloSort direcory (ksDir)
% should be done AFTER manula curation stage in Phy

%% add necessary paths to toolboxes:

paths = addPathsForSpikeSorting;

%% list of paths to raw ephys files
ksDirList = {...
     'Y:\LAB PROJECTS\scInactivateAndRecord\data\ramsey20191126\';...
% 'Y:\LAB PROJECTS\scDualVprobe\data\ram20190805a\'; ...
% 'Y:\Leor\170421\'; ...
% 'Y:\Leor\170510\'; ...
% 'Y:\Leor\170522\'; ...
    }; 

nFiles = numel(ksDirList);

%% (1) load sp, (2) mkfigs, (3) convert to su, (4) mkfigs

for iF = 1:nFiles
    tic
    ksDir = fullfile(ksDirList{iF}, 'kiloSorted');
    
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(['~~~~~  ' ksDir])
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
   % load sp:
    sp = getSp(ksDir);
   
   % make sure there's a figures folder:
   if ~exist(fullfile(ksDir, 'figures'), 'dir')
       mkdir(fullfile(ksDir, 'figures'));
   end
   
   % set options:
   opts.saveFigs = 1;
   opts.dirFigs = ksDir;
   
   % make dem figures and save:
   mkfig.waveform_overChannels(sp, opts);
   mkfig.waveformOverChannels_perCluster(sp, opts);
   mkfig.waveformAndSpikeCount_overChannels(sp, opts)
   su = sp2su(sp, ksDir);
   mkfig.unitSummary(su, opts);
   
    close all
   timePassed(iF) = toc;
end


%%


