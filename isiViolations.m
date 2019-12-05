function isiV = isiViolations(IDs, varargin)

if nargin > 1
    resultsDirectory = varargin{1};
else
    resultsDirectory = uigetdir();
end

%   [cgs, uQ, cR, isiV] = computeAllMeasures(resultsDirectory)
%
% katz comments:
% INPUT:
%   resultsDirectory - directory with kilo output
% OUTPUT:
%   cgs  - cluster groups (0=noise; 1=MUA; 2=Good; 3=Unsorted)
%   uQ   - unitQuality AKA isolation distance
%   cR   - contamination rate = the proportion of spikes inside the cluster 
%          boundary that aren't from the cluster (false positive rate)
%   isiV - isi Violations = the estimated false positive rate of your spike
%          train, based on the rate of refractory period violations.

% clusterPath = fullfile(resultsDirectory, 'cluster_groups.csv');
% spikeClustersPath = fullfile(resultsDirectory,'spike_clusters.npy');
% spikeTemplatesPath = fullfile(resultsDirectory,'spike_templates.npy');

% if exist(clusterPath, 'file')
%     % ie a csv with sorted spike info exists
%     [cids, cgs] = readClusterGroupsCSV(clusterPath);
% elseif exist(spikeClustersPath, 'file')
%     % ie a csv doesn't exist, but we can get info from the NPYs:
%     clu = readNPY(spikeClustersPath);
%     cgs = 3*ones(size(unique(clu))); % set all to "unsorted"
% else
%     clu = readNPY(spikeTemplatesPath);
%     cgs = 3*ones(size(unique(clu))); % all unsorted
% end

% refractory duration & minimum isi definition
refDur = 0.0015;
minISI = 0.0005;

nIDs = length(IDs);

fprintf(1, 'loading data...\n');
%% Precompute the locationsn of files to be loaded
pcFeaturesPath      = fullfile(resultsDirectory,'pc_features.npy');
% pcFeaturesIndPath   = fullfile(resultsDirectory,'pc_feature_ind.npy');
spikeClustersPath   = fullfile(resultsDirectory,'spike_clusters.npy');
% spikeTemplatesPath  = fullfile(resultsDirectory,'spike_templates.npy');
spikeTimesPath      = fullfile(resultsDirectory,'spike_times.npy');
paramsPath          = fullfile(resultsDirectory,'params.py');

%% Main code.
try
%     pc_features = readNPY(pcFeaturesPath); % features of each spike
catch me
    if ~exist(pcFeaturesPath, 'file')
        fprintf(1, 'PC Features loading failed. File does not exist.\n');
    else
        fprintf(1, 'PC Features loading failed. You may need to clone the npy-matlab repo and add to path.\n');
    end
    rethrow me
end
% pc_feature_ind = readNPY(pcFeaturesIndPath); % the (small) subset of channels corresponding to each template
% pc_feature_ind = pc_feature_ind + 1;   % because in Python indexes start at 0

if exist(spikeClustersPath,'file')
%     fprintf('building features matrix from clusters/templates\n')
    spike_clusters = readNPY(spikeClustersPath);
%     spike_clusters = spike_clusters; % because in Python indexes start at 0
        
    % now we have to construct a new pc_features that has the features
    % arranged according to *cluster* rather than template. 
%     spike_templates = readNPY(spikeTemplatesPath);
%     spike_templates = spike_templates+1; % !!!!!!!!!!!!!!!!!!! ??????
    
    clusterIDs = IDs;
%     nClusters = length(clusterIDs);
%     nSpikes = length(spike_clusters);
%     nFet = 4; nFetPerChan = size(pc_features,2);
%     nTemplates = size(pc_feature_ind,1);
%     
%     newFet = zeros(nSpikes, nFetPerChan, nFet);
%     newFetInds = zeros(nClusters, nFet);
    %tempNums = 1:nTemplates;
    
    spike_times = readNPY(spikeTimesPath);
params = readKSparams(paramsPath);
spike_times = double(spike_times)/params.sample_rate;
% 
fprintf(1, 'computing ISI violations\n');
% 
    for c = 1:nIDs
%         fprintf(1, '%d/%d\n', c, length(clusterIDs))
%         thisID = clusterIDs(c);
%         
%         theseSpikes = spike_clusters==thisID;
%         theseTemplates = spike_templates(theseSpikes);
%         
%         [inclTemps, inst] = countUnique(theseTemplates);
%         
%         thisTemplate = inclTemps(inst==max(inst),1);
%         
%         theseChans = pc_feature_ind(thisTemplate,1:nFet);
%         
%         
%         
%         newFetInds(c,:) = theseChans;
        
        %subPCFetInd = pc_features(theseSpikes,:,:);
        
                
%         for f = 1:nFet
%             thisChanInds = pc_feature_ind==theseChans(f);
%             [chanInds,tempsWithThisChan] = find(thisChanInds');
%             %spikesWithThisFet = ismember(theseTemplates, tempsWithThisChan);
%                         
%             inclTempsWithThisFet = find(ismember(inclTemps, tempsWithThisChan));
%             for t = 1:numel(inclTempsWithThisFet)
%                 thisSubTemp = inclTemps(inclTempsWithThisFet(t));
%                 thisTfetInd = chanInds(tempsWithThisChan==thisSubTemp);
%                 newFet(theseSpikes&spike_templates==thisSubTemp,:,f) = ...
%                     pc_features(theseSpikes&spike_templates==thisSubTemp,:,thisTfetInd);
%             end
%             
%             
%         end
        
        [fpRate, numViolations] = ISIViolations(spike_times(spike_clusters==clusterIDs(c)), minISI, refDur);
        isiV(c) = fpRate;
        nSpikes = sum(spike_clusters==clusterIDs(c));
        fprintf(1, 'cluster %3d: %d viol (%d spikes), %.2f estimated FP rate\n', clusterIDs(c), numViolations, nSpikes, fpRate);
        
    end
    
%     pc_features = newFet;
%     all_pc_feature_ind = pc_feature_ind;
%     pc_feature_ind = newFetInds;
else
    fprintf(1, 'warning, spike_clusters does not exist, using spike_templates instead\n');
%     spike_clusters = readNPY(spikeTemplatesPath); % template # corresponding to each spike
%     spike_clusters = spike_clusters +1; % because in Python indexes start at 0
end



% assert(numel(size(pc_features)) == 3)

% fprintf(1, 'computing cluster qualities...\n');
% 
% % - spike_clusters is 1 x nSpikes
% % - pc_features is nSpikes x nPCsPerChan x nInclChans
% % - pc_feature_ind is nClusters x nInclChans (sorted in descending order of
% % relevance for this template)
% % - fetN is an integer, the number of features to
% 
% if ~exist('fetNchans', 'var')
%     fetNchans = min(4, size(pc_feature_ind,2)); % number of channels to use
% end
% nFetPerChan = size(pc_features,2);
% fetN = fetNchans*nFetPerChan; % now number of features total
% % pc_features = reshape(pc_features, size(pc_features,1), []);
% 
% N = numel(spike_clusters);
% assert(fetNchans <= size(pc_features, 3) && size(pc_features, 1) == N , 'bad input(s)')
% 
% unitQuality = zeros(size(clusterIDs));
% contaminationRate = zeros(size(clusterIDs));
% 
% allClusters = unique(spike_clusters);
% nAllClusters = length(allClusters);
% 
% fprintf('%12s\tQuality\tContamination\n', 'ID'); % comment to suppress printing out the intermediate results
% for c = 1:nIDs
%     
%     theseSp = spike_clusters == clusterIDs(c);
%     n = sum(theseSp); % #spikes in this cluster
%     if n < fetN || n >= N/2
%         % cannot compute mahalanobis distance if less data points than
%         % dimensions or if > 50% of all spikes are in this cluster
%         unitQuality(c) = 0;
%         contaminationRate(c) = NaN;
%         continue
%     end
%     
%     fetThisCluster = reshape(pc_features(theseSp,:,1:fetNchans), n, []);
%     
%     % now we need to find other spikes that exist on the same channels
%     theseChans = pc_feature_ind(c,1:fetNchans);
% 
%     % for each other cluster, determine whether it has at least one of
%     % those channels. If so, add its spikes, with its features put into the
%     % correct places
%     nInd = 1; fetOtherClusters = zeros(0,size(pc_features,2),fetNchans);
%     for c2 = 1:nAllClusters
%         if c2~=c
%             chansC2Has = all_pc_feature_ind(c2,:);
%             for f = 1:length(theseChans)
%                 
%                 if ismember(theseChans(f), chansC2Has)
%                     
%                     theseOtherSpikes = spike_clusters==allClusters(c2);
%                     thisCfetInd = find(chansC2Has==theseChans(f),1);
%                     fetOtherClusters(nInd:nInd+sum(theseOtherSpikes)-1,:,f) = ...
%                         pc_features(theseOtherSpikes,:,thisCfetInd);
%                 end                
%                 
%             end
%             if any(ismember(chansC2Has, theseChans))
%                 nInd = nInd+sum(theseOtherSpikes);
%             end
%         end
%     end
%     
%                 
%                 
%                 
%     
%     fetOtherClusters = reshape(fetOtherClusters, size(fetOtherClusters,1), []);
%     
%     [uQ, cR] = maskedClusterQualityCore(fetThisCluster, fetOtherClusters);
%     
%     unitQuality(c) = uQ;
%     contaminationRate(c) = cR;
%     
%     fprintf('cluster %3d: \t%6.1f\t%6.2f\n', clusterIDs(c), unitQuality(c), contaminationRate(c)); % comment to suppress printing out the intermediate results
%     
%     if uQ>1000
%         keyboard;
%     end
    



%%
% %% Precompute the locationsn of files to be loaded
% spikeClustersPath = fullfile(resultsDirectory,'spike_clusters.npy');
% spikeTemplatesPath = fullfile(resultsDirectory,'spike_templates.npy');
% spikeTimesPath= fullfile(resultsDirectory,'spike_times.npy');
% paramsPath= fullfile(resultsDirectory,'params.py');
% 
% %% 
% 

% 
% fprintf(1, 'loading data for ISI computation\n');
% if exist(spikeClustersPath)
%     spike_clusters = readNPY(spikeClustersPath);
% else
%     spike_clusters = readNPY(spikeTemplatesPath);
% end
% spike_clusters = spike_clusters + 1; % because in Python indexes start at 0
% 
% spike_times = readNPY(spikeTimesPath);
% params = readKSparams(paramsPath);
% spike_times = double(spike_times)/params.sample_rate;
% 
% fprintf(1, 'computing ISI violations\n');
% 
% clusterIDs = IDs;
% isiV = zeros(1,numel(clusterIDs));
% for c = 1:nIDs
%     
%     [fpRate, numViolations] = ISIViolations(spike_times(spike_clusters==clusterIDs(c)), minISI, refDur);
%     isiV(c) = fpRate;
%     nSpikes = sum(spike_clusters==clusterIDs(c));    
%     fprintf(1, 'cluster %3d: %d viol (%d spikes), %.2f estimated FP rate\n', clusterIDs(c), numViolations, nSpikes, fpRate);
%     
% end

end


% %%
% if iPlot==1
%         title('info')
%     end
%     hold on
%     axis([0 1 0 1])
%     yPosStart = 0.1;
%     if isfield(thisUnit, 'times')
%         text(0.1, yPosStart, ['nSpikes = ', sprintf('%.0f', numel(thisUnit.times))]); 
%         yPosStart = yPosStart+0.2;
%     end
%     if isfield(thisUnit, 'uQ')
%         text(0.1, yPosStart, ['uQ = ', sprintf('%.2f', thisUnit.uQ)]); 
%         yPosStart = yPosStart+0.2;
%     end
%     if isfield(thisUnit, 'cR')
%         text(0.1, 0.7, ['cR = ', sprintf('%.3f', thisUnit.cR)]);
%         yPosStart = yPosStart+0.2;
%     end
%     if isfield(thisUnit, 'isiV')
%         text(0.1, 0.5, ['isiV = ', sprintf('%.3f', thisUnit.isiV)]);
%         yPosStart = yPosStart+0.2;
%     end
%     if isfield(thisUnit, 'isolation')
%         text(0.1, 0.5, ['isolation = ', sprintf('%.2f', thisUnit.isolation)]);
%         yPosStart = yPosStart+0.2;
%     end
%     
