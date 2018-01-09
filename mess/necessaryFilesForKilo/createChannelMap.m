% dirDataset = '~/Dropbox/_transfer_big_data/20170615/data/convertedToRaw';
%  create a channel map file

Nchannels   = 24;
connected   = true(Nchannels, 1);
chanMap     = 1:Nchannels;
chanMap0ind = chanMap - 1;

% sterotrode geometry
% xSpacing = 200; % phy visualization is weird with x axis so I'm giving it 200 (despite real value being 50).
% ySpacing = 50;
% xcoords = xSpacing * repmat([0 1]', Nchannels/2,1);
% ycoords = ySpacing * kron(((Nchannels/2):-1:1)-1, [1 1])'; % using kron for stereo geometry
% linear geometry:
xcoords   = ones(Nchannels,1);
ycoords   = [1:Nchannels]';
kcoords = ones(Nchannels,1); 

fs = 44000; % sampling frequency
save(fullfile(dirDataset, 'chanMap.mat'), ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')

%%

% kcoords is used to forcefully restrict templates to channels in the same
% channel group. An option can be set in the master_file to allow a fraction 
% of all templates to span more channel groups, so that they can capture shared 
% noise across all channels. This option is

% ops.criterionNoiseChannels = 0.2; 

% if this number is less than 1, it will be treated as a fraction of the total number of clusters

% if this number is larger than 1, it will be treated as the "effective
% number" of channel groups at which to set the threshold. So if a template
% occupies more than this many channel groups, it will not be restricted to
% a single channel group. 