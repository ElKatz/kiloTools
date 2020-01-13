function sp = getSp(ksDir, varargin)
%   sp = getSp(ksDir, varargin)
%
% get a spike struct from a kiloSort output directory (ksDir)
% The function reads npy files from the ksDir to populate the sp struct
% with all you ever dreamt of and more.
% 
% It gets all spike info and computes certain spike metrics such as ISI
% violation and cluster quality using routines stolen from CortexLab, using
% the sqKilosort.computeAllMeasures function.
%
% INPUT:
%   ksDir - path to directory of kiloSort output
%   varargin:
%   (see default wihin code)
%       waves   - if true, adds individual waveforms to the strcut (time
%                 consuming, so default is false)
%       exNoise - if true, excludes noise clusters (as defined in phy)
%       exMu    - if true, excludes multiunit clusters (as defined in phy)
%       loadPCs - if true, loads PCs.
%       waveWinT - [2x1] vector of time before & after spike time to extract
%                 the full waveform. in ms.
%
% OUTPUT:
%   sp - struct with all the goodies you need to access your spike data. 
% 

%% info:

%   dat_path:         your filename e.g. 'expt52.dat'
%   n_channels_dat:   e.g. 24
%   dtype:            e.g. 'int16'
%   offset:           offset between recording start time and sort time.
%   sample_rate:      sample rate of your DAQ system, e.g. 40000
%   hp_filtered:      0
%   dat_folder:       folder for your filename e.g. '/Users/leorkatz/Dropbox/Code/spike_sorting/pilot_dat?'
%   ss:             sample number at which there was a spike from any cluster 
%   st:             [nTotalSpikes x 1] !!!! 2do !!! (* naming: st_s to indicate that it
%       is in seconds.
%   spikeTemplates: [nTotalSpikes x 1] the template number associated with
%       each spike.
%   clu: [nTotalSpikes x 1] the cluster number associated with
%       each spike.
%   tempScalingAmps: [nTotalSpikes x 1] the template scaling amplitude for
%   each spike
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
% medWfs: [nClusters × nCh × nWaveformSamples] median waveform on channel
% which had the largest waveform amplitude. 
%     




%% defaults: 

p = inputParser;
p.addOptional('waves', false);
p.addOptional('exNoise', true);
p.addOptional('exMu', true);
p.addOptional('loadPCs', false);
p.addOptional('waveWinT', [-300 900]);
p.addOptional('medWave', true)
p.addOptional('visualize', false)
p.addOptional('save', true)
p.parse(varargin{:});

%% Store the info:

% load thee info struct
load(fullfile(ksDir, 'convertInfo.mat'));
info.ksDir          = ksDir;
info.depth          = [];
info.dTip2lowestCh  = [];
info.b              = [];

%% LOAD UP DATA FROM npy & csv FILES:

% load the sampsToSecMap:
load(fullfile(ksDir, 'sampsToSecsMap.mat'));

% load spike data from npy:
sp              = loadParamsPy(fullfile(ksDir, 'params.py'));
sp.info         = info;
spikeTimesSamps = readNPY(fullfile(ksDir, 'spike_times.npy'));
spikeTimesSecs  = sampsToSecsMap(spikeTimesSamps);

%% convert samples tp seconds:
% (AKA convert spikeTimesSamps to spikeTimeSecs):

% 
% % read the strobed word info (values & time stamps).
% out.eventInfo = PL2EventTs(out.PL2fileName, 'Strobed');
% out.RSTARTInfo = PL2EventTs(out.PL2fileName, 'RSTART');
% out.RSTOPInfo = PL2EventTs(out.PL2fileName, 'RSTOP');
% 
% % must read in a spike channel to construct the "timestamp map" from
% % samples (kilosort) to time in seconds.
% ad = PL2Ad(out.PL2fileName, 'SPKC01');
% 
% % place to store the "map" from samples to seconds.
% sampsToSecsMap = zeros(sum(ad.FragCounts),1);
% 
% % sample duration
% sampDur = 1/ad.ADFreq;
% 
% % how many fragments of recording?
% nFrags = length(ad.FragTs);
% currentSample = 1;
% for i = 1:nFrags
%     chunkIndex = currentSample:(currentSample + ad.FragCounts(i) - 1);
%     timeStamps = ad.FragTs(i) + (0:(ad.FragCounts(i)-1))*sampDur;
%     sampsToSecsMap(chunkIndex) = timeStamps;
%     currentSample = chunkIndex(end)+1;
% end
% 
% % map samples to timeStamps
% out.spikeTimesSecs = sampsToSecsMap(out.spikeTimesSamps);
% spikeTimesSecs = double(spikeTimesSamps)/sp.sample_rate;


%%



spikeTemplates  = readNPY(fullfile(ksDir, 'spike_templates.npy')); % note: zero-indexed
if exist(fullfile(ksDir, 'spike_clusters.npy'), 'file')
    spikeClusters         = readNPY(fullfile(ksDir, 'spike_clusters.npy'));
else
    spikeClusters         = spikeTemplates;
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
    csvFile         = fullfile(ksDir, 'cluster_groups.csv');
    [clusterId, clusterScore]  = readClusterGroupsCSV(csvFile);
else
    error('why no csv file? you baaad')
end

% if you wish to exlcude noise or multiunit, this is where it happens:
cIdExclude = [];
if p.Results.exNoise
    cIdExclude = [cIdExclude, clusterId(clusterScore==0)];
end
if p.Results.exMu
    cIdExclude = [cIdExclude, clusterId(clusterScore==1)];
end

% exclude'em:
spikeTimesSamps = spikeTimesSamps(~ismember(spikeClusters, cIdExclude));
spikeTimesSecs  = spikeTimesSecs(~ismember(spikeClusters, cIdExclude));
spikeTemplates  = spikeTemplates(~ismember(spikeClusters, cIdExclude));
tempScalingAmps = tempScalingAmps(~ismember(spikeClusters, cIdExclude));
spikeClusters   = spikeClusters(~ismember(spikeClusters, cIdExclude));
clusterScore    = clusterScore(~ismember(clusterId, cIdExclude));
clusterId       = clusterId(~ismember(clusterId, cIdExclude));

if p.Results.loadPCs
    pcFeat = pcFeat(~ismember(spikeClusters, cIdExclude), :,:);
    %pcFeatInd = pcFeatInd(~ismember(cids, cidsExclude),:);
end


% get the last pieces of data:
coords  = readNPY(fullfile(ksDir, 'channel_positions.npy'));
xcoords = coords(:,1);
ycoords = coords(:,2);
temps   = readNPY(fullfile(ksDir, 'templates.npy'));
winv    = readNPY(fullfile(ksDir, 'whitening_mat_inv.npy'));

% and pack it up:
sp.spikeTimesSamps  = spikeTimesSamps;
sp.spikeTimesSecs   = spikeTimesSecs;
sp.spikeTemplates   = spikeTemplates;
sp.spikeClusters    = spikeClusters;
sp.tempScalingAmps  = tempScalingAmps;
sp.clusterScore     = clusterScore;
sp.clusterId        = clusterId;
sp.clusterStr       = arrayfun(@(x) sprintf('%0.3d', x), clusterId, 'UniformOutput', 0); % for easy figure titles
sp.nClusters        = numel(sp.clusterId);
sp.nChannels        = numel(xcoords);
sp.xcoords          = xcoords;
sp.ycoords          = ycoords;
sp.temps            = temps;
sp.winv             = winv;
sp.pcFeat           = pcFeat;
sp.pcFeatInd        = pcFeatInd;
sp.wf               = [];

%% Sort quality meausres:
disp('Computing sort quality measures...')

%   uQ   - unitQuality AKA isolation distance
%   cR   - contamination rate = the proportion of spikes inside the cluster 
%          boundary that aren't from the cluster (false positive rate)
%   isiV - isi Violations = the estimated false positive rate of your spike
%          train, based on the rate of refractory period violations.

[~, sp.uQ, sp.cR] = sqKilosort.maskedClusterQuality(ksDir);
[sp.isiV_fpRate, sp.isiV_rate] = compute_isiViolations(ksDir);
   

%% waves:
% if you wish to get waveforms (time consuming), this is where it happens:
if p.Results.waves
    disp('Retreiving all (yes all) waveforms. This might take a while...')
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
    spikeTimesSamps      = double(sp.spikeTimesSamps);
    %
    % get raw data around spiketimes (in samples) and populate 'wv' of size
    % [nCh, nWaveSamples, nSpikes]
    wf = zeros(nCh, numel(win), numel(spikeTimesSamps), 'int16');
    % for each spikes:
    for ii = 1:length(spikeTimesSamps)
        spkwin = spikeTimesSamps(ii) + win;
        wf(:,:,ii) = dat(:,spkwin);
    end
    sp.wv = wf;
end

%% median waveform
% kinda just hacking this section together so probably not the most
% efficient or clean...
if p.Results.medWave
    disp('Retreiving median waveforms...')
    load(fullfile(ksDir, 'ops.mat'));
%     d = dir(ops.fbinary);
    d = dir([info.ksDir '/*.dat']);
    nSamp = d.bytes/2/sp.n_channels_dat;
    dataSize = [sp.n_channels_dat nSamp];
    chanMap = readNPY(fullfile(ksDir, 'channel_map.npy'));
    gain = 0.6/512/500*1e6; % raw file units to uV ***SHOULD BE RIG SPECIFIC. NEED TO DO THIS...
    
    % median wf per cluster per channel, size: [nClu, nCh, nSamps]
%     sp.medWfs = extractMedianWFs(sp.spikeClusters, sp.spikeTimesSecs, sp.sample_rate, fullfile(ksDir, sp.dat_path), sp.dtype, dataSize, chanMap, gain);
    sp.medWfs = extractMedianWFs(sp.spikeClusters, sp.spikeTimesSamps, sp.sample_rate, fullfile(ksDir, sp.dat_path), sp.dtype, dataSize, chanMap, gain, 'samples');
    % median wf per cluster on peak amplitude channel, size: [nClu, nSamps]
    for iS = 1:sp.nClusters 
        medWfPerCh              = squeeze(sp.medWfs(iS,:,:))';
        medWfAmpPerCh           = max(medWfPerCh) - min(medWfPerCh);
        sp.peakCh(iS)           = find(medWfAmpPerCh == max(medWfAmpPerCh), 1);
        sp.medWfOnPeakCh(iS,:)  = medWfPerCh(:, sp.medWfPeakCh(iS));
    end
    
    if p.Results.visualize
        mkfig.medWfPerChannel(sp)
    end
end

%% save sp

if p.Results.save
    save(fullfile(ksDir, 'sp.mat'), '-struct', 'sp')
    disp('Done saving ''sp''')
else
    disp('not saving')
end


disp('-------------------')
disp('DONE! Enjoy your sp')
disp('-------------------')




