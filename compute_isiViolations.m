function [isiV_fpRate, isiV_rate] = compute_isiViolations(resultsDirectory)

%% Precompute the locationsn of files to be loaded
spikeClustersPath = fullfile(resultsDirectory,'spike_clusters.npy');
spikeTemplatesPath = fullfile(resultsDirectory,'spike_templates.npy');
spikeTimesPath= fullfile(resultsDirectory,'spike_times.npy');
paramsPath= fullfile(resultsDirectory,'params.py');

%% 

refDur = 0.001;     % estiamtion of refractory period duration (lnk: was 0.0015)
minISI = 0.0005;    % min possible ISI given waveform length (lnk: this is only used to compute the false positive rate);

fprintf(1, 'loading data for ISI computation\n');
if exist(spikeClustersPath)
    spike_clusters = readNPY(spikeClustersPath);
else
    spike_clusters = readNPY(spikeTemplatesPath);
end
spike_clusters = spike_clusters + 1; % because in Python indexes start at 0

spike_times = readNPY(spikeTimesPath);
params = readKSparams(paramsPath);
spike_times = double(spike_times)/params.sample_rate;

fprintf(1, 'computing ISI violations\n');

clusterIDs  = unique(spike_clusters);
isiV_fpRate = zeros(numel(clusterIDs),1);
isiV_rate   = zeros(numel(clusterIDs),1);

for c = 1:numel(clusterIDs)
    
    [fpRate, numViolations] = ISIViolations(spike_times(spike_clusters==clusterIDs(c)), minISI, refDur);
    isiV_fpRate(c) = fpRate;
    nSpikes = sum(spike_clusters==clusterIDs(c));    
    isiV_rate(c) = numViolations/nSpikes;
    fprintf(1, 'cluster %3d: %d violations out of %d spikes. %.2f isiV_rate, %.2f isiV_fpRate\n', ...
        clusterIDs(c), numViolations, nSpikes, isiV_rate(c), isiV_fpRate(c));
    
end
