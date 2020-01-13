function [paths] = addPathsForSpikeSorting
%
% Different machines have different paths... This script adds the paths for
% your spike sorting needs given a particular machine.
% Fee free to add your own machine to the list...
%
% Necessary paths:
%   kiloTools
%   KiloSort-master
%   npy-matlab
%
% INPUT
%   none
% OUTPUT
%   paths - struct of folder locations

%% Define paths

% get hostname:
[~, hostName] = system('hostname');

% addpaths
if contains(hostName, 'DESKTOP-KEJGC64', 'IgnoreCase', 1)
    % Dell spike sorter in NIH krauzlis rig space 
    paths.lnkToolbox    = 'D:\Code\Toolboxes\0-lnkToolbox';
    paths.kiloTools     = 'D:\Code\Toolboxes\kiloTools';
    paths.kiloSort      = 'D:\Code\Toolboxes\KiloSort-master';
    paths.npymatlab     = 'D:\Code\Toolboxes\npy-matlab';
    
elseif contains(hostName, 'LA-CPS828317MN-Huk-2.local', 'IgnoreCase', 1) || ... % Leor UT MBP
        contains(hostName, 'NEIK2A79LK07A', 'IgnoreCase', 1) || ... % Leor NIH iMac
        contains(hostName, 'NEIK2A50LK49B', 'IgnoreCase', 1) % Leor NIH MBP
    
    paths.lnkToolbox    = '~/Dropbox/Code/spike_sorting/0-lnkToolbox';
    paths.kiloTools      = '~/Dropbox/Code/spike_sorting/toolboxes/kiloTools';
    paths.spikes         = '~/Dropbox/Code/spike_sorting/toolboxes/spikes';
    paths.sortingQuality = '~/Dropbox/Code/spike_sorting/toolboxes/sortingQuality';
    paths.kiloSort       = '~/Dropbox/Code/spike_sorting/packages/KiloSort';
    paths.npymatlab      = '~/Dropbox/Code/spike_sorting/toolboxes/npy-matlab';
    paths.plexonSdk      = '~/Dropbox/Code/Tools/Plexon Offline SDKs/Matlab Offline Files SDK';
    
else
    error('Unrecognized hostname. Could not add necessary paths for sorting')
end

%% add paths
flds = fieldnames(paths);
for iF = 1:numel(flds)
    addpath(genpath(paths.(flds{iF})));
    disp(['added - ' flds{iF}]);
end

