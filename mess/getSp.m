function sp = getSp(ksDir, varargin)
%   sp = getSp(ksDir, varargin)
%
% get a spike struct from a kiloSort output directory (ksDir)
% The function reads npy files from the ksDir to populate the sp struct
% with all you ever dreamt of and more.
% INPUT:
%   ksDir - path to directory of kiloSort output
%   optional:
%       waves   - if true, adds individual waveforms to the strcut (time
%                 consuming, so default is false)
%       exNoise - if true, excludes noise clusters (as defined in phy)
%       exMu    - if true, excludes multiunit clusters (as defined in phy)
%       loadPCs - if true, loads PCs.
%       waveWinT - [2x1] vector of time before & after spike time to extract
%                 the full waveform. in ms.


%% info:
% dat_path: your filename e.g. 'expt52.dat'
% n_channels_dat: e.g. 24
% dtype: e.g. 'int16'
% offset: offset between recording start time and sort time.
% sample_rate: sample rate of your DAQ system, e.g. 40000
% hp_filtered: 0
% dat_folder: folder for your filename e.g. '/Users/leorkatz/Dropbox/Code/spike_sorting/pilot_dat?'
%%% the following are in size [nTotalSpikes x 1] i.e. for each spike
% ss:  sample number at which there was a spike from any cluster 
% st: [nTotalSpikes x 1] !!!! 2do !!! (* naming: st_s to indicate that it
% is in seconds.
% spikeTemplates: [nTotalSpikes x 1] the template number associated with
% each spike.
% clu: [nTotalSpikes x 1] the cluster number associated with
% each spike.
% tempScalingAmps: [nTotalSpikes x 1] the template scaling amplitude for
% each spike
%%% the following are in size [nClusters x 1] i.e. for each sorted cluster
% cgs: [nClusters x 1] group number of each cluster e.g. [2 1 1 0 2] where 0=noise,
% 1=multiunit, 2=good.
% cids: [nClusters x 1] cluster id for each cluster e.g. [1 7 22 23 38]
% nClu: number of clusters
% nCh: number of channels
% xcoords: [nCh×1] x locations of recording channels
% ycoords: [nCh×1] y locations of recording channels
% temps: [nTemplates × nTemplateSamples × nCh] template shapes used for
% kilosort
% winv: [24×24 single]  !!!!! 2do 
% pcFeat: [] !!!!! 2do 
% pcFeatInd: [] !!!!! 2do 
% wv: [nCh × nWaveformSamples × nTotalSpikes] full waveform within a window
% (window size defined in 'getSp.m')
% medWFs: [nClusters × nCh × nWaveformSamples] median waveform on channel
% which had the largest waveform amplitude. 
%     




%%
p = inputParser;
p.addOptional('waves', false);
p.addOptional('exNoise', true);
p.addOptional('exMu', true);
p.addOptional('loadPCs', false);
p.addOptional('waveWinT', [-300 900]);
p.addOptional('medWave', true)
p.addOptional('visualize', false)
p.parse(varargin{:});


%% LOAD UP DATA FROM npy & csv FILES:

% load spike data from npy:
sp              = loadParamsPy(fullfile(ksDir, 'params.py'));
sp.dat_folder   = ksDir;
ss              = readNPY(fullfile(ksDir, 'spike_times.npy'));
st              = double(ss)/sp.sample_rate;
spikeTemplates  = readNPY(fullfile(ksDir, 'spike_templates.npy')); % note: zero-indexed
if exist(fullfile(ksDir, 'spike_clusters.npy'), 'file')
clu         = readNPY(fullfile(ksDir, 'spike_clusters.npy'));
else
    clu         = spikeTemplates;
end
tempScalingAmps = readNPY(fullfile(ksDir, 'amplitudes.npy'));
if p.Results.loadPCs
    pcFeat      = readNPY(fullfile(ksDir,'pc_features.npy')); % nSpikes x nFeatures x nLocalChannels
    pcFeatInd   = readNPY(fullfile(ksDir,'pc_feature_ind.npy')); % nTemplates x nLocalChannels
else
    pcFeat      = [];
    pcFeatInd   = [];
end

% get the phy output data from csv:
if exist(fullfile(ksDir, 'cluster_groups.csv'), 'file')
    cgsFile     = fullfile(ksDir, 'cluster_groups.csv');
    [cids, cgs] = readClusterGroupsCSV(cgsFile);
else
    error('why no csv file? you baaad')
end

% if you wish to exlcude noise or multiunit, this is where it happens:
cidsExclude = [];
if p.Results.exNoise
    cidsExclude = [cidsExclude, cids(cgs==0)];
end
if p.Results.exMu
    cidsExclude = [cidsExclude, cids(cgs==1)];
end

% exclude'em:
ss              = ss(~ismember(clu, cidsExclude));
st              = st(~ismember(clu, cidsExclude));
spikeTemplates  = spikeTemplates(~ismember(clu, cidsExclude));
tempScalingAmps = tempScalingAmps(~ismember(clu, cidsExclude));
clu             = clu(~ismember(clu, cidsExclude));
cgs             = cgs(~ismember(cids, cidsExclude));
cids            = cids(~ismember(cids, cidsExclude));

if p.Results.loadPCs
    pcFeat = pcFeat(~ismember(clu, cidsExclude), :,:);
    %pcFeatInd = pcFeatInd(~ismember(cids, cidsExclude),:);
end


% get the last pieces of data:
coords  = readNPY(fullfile(ksDir, 'channel_positions.npy'));
xcoords = coords(:,1);
ycoords = coords(:,2);
temps   = readNPY(fullfile(ksDir, 'templates.npy'));
winv    = readNPY(fullfile(ksDir, 'whitening_mat_inv.npy'));

% and pack it up:
sp.ss               = ss;
sp.st               = st;
sp.spikeTemplates   = spikeTemplates;
sp.clu              = clu;
sp.tempScalingAmps  = tempScalingAmps;
sp.cgs              = cgs;
sp.cids             = cids;
sp.nClu             = numel(sp.cids);
sp.nCh              = numel(xcoords);
sp.xcoords          = xcoords;
sp.ycoords          = ycoords;
sp.temps            = temps;
sp.winv             = winv;
sp.pcFeat           = pcFeat;
sp.pcFeatInd        = pcFeatInd;
sp.wv               = [];

%% waves:
% if you wish to get waveforms (time consuming), this is where it happens:
if p.Results.waves
    
    % load in raw data
    load(fullfile(ksDir, 'ops.mat'));
    fid     = fopen(ops.fbinary, 'r');
    chMap   = load(fullfile(ksDir, 'chanMap.mat'));
    nCh     = numel(chMap.chanMap);
    dat     = fread(fid, [nCh inf], '*int16');
    dat     = dat(chMap.connected,:); % only take "connected" channels
    fclose(fid);
    
    
    % define window within which to take the waveform:
    winT    = p.Results.waveWinT .* 1e-6; % convert to seconds
    winS    = winT .* sp.sample_rate; % in samples
    win     = floor(winS(1):winS(end));
    
    % spike_times.npy is actually spike sample (ss) of every spike. (to
    % convert to time simply divide by Fs)
    %     ss      = readNPY(fullfile(ksDir,  'spike_times.npy'));
    ss      = double(sp.ss);
    %
    % get raw data around spiketimes (in samples) and populate 'wv' of size
    % [nCh, nWaveSamples, nSpikes]
    wv = zeros(nCh, numel(win), numel(ss), 'int16');
    % for each spikes:
    for ii = 1:length(ss)
        spkwin = ss(ii) + win;
        wv(:,:,ii) = dat(:,spkwin);
    end
    sp.wv = wv;
end

%% median waveform
% kinda just hacking this section together so probably not the most
% efficient or clean...
if p.Results.medWave
    load(fullfile(ksDir, 'ops.mat'));
    d = dir(ops.fbinary);
    nSamp = d.bytes/2/sp.n_channels_dat;
    dataSize = [sp.n_channels_dat nSamp];
    chanMap = readNPY(fullfile(ksDir, 'channel_map.npy'));
    gain = 0.6/512/500*1e6; % raw file units to uV ***SHOULD BE RIG SPECIFIC. NEED TO DO THIS...
    sp.medWFs = extractMedianWFs(sp.clu, sp.st, sp.sample_rate, fullfile(ksDir,sp.dat_path), sp.dtype, dataSize, chanMap, gain);
    if p.Results.visualize
        mkfig_medWF(sp)
    end
end





