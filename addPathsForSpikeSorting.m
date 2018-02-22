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
if contains(hostName, 'lsr-rjk-mata', 'IgnoreCase', 1)
    % Dell spike sorter in Leor office
    paths.kiloTools = 'C:\EPHYS\Code\Toolboxes\kiloTools';
    paths.kiloSort  = 'C:\EPHYS\Code\Toolboxes\KiloSort-master';
    paths.npymatlab = 'C:\EPHYS\Code\Toolboxes\npy-matlab';
    
elseif contains(hostName, 'LA-CPS828317MN-Huk-2.local', 'IgnoreCase', 1)
    % Leor MacBookPro
    paths.kiloTools = '~/Dropbox/Code/spike_sorting/toolboxes/kiloTools';
    paths.kiloSort  = '~/Dropbox/Code/spike_sorting/packages/KiloSort';
    paths.npymatlab = '~/Dropbox/Code/spike_sorting/toolboxes/npy-matlab';
    
else
    error('Unrecognized hostname. Could not add necessary paths for sorting')
end

%% add paths
addpath(genpath(paths.kiloTools))
addpath(genpath(paths.kiloSort))
addpath(genpath(paths.npymatlab))
disp('Paths for spikesorting-- added!')