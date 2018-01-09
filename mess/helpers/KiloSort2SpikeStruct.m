function out = KiloSort2SpikeStruct(rez)

% load in raw data
fid = fopen(rez.ops.fbinary, 'r');
NchanTOT = rez.ops.NchanTOT;
dat = fread(fid, [NchanTOT inf], '*int16');
fclose(fid);

% organize data with chanMap, remove unconnected channels
dat = dat(rez.ops.chanMap(rez.connected),:);
win = [-50:50];

% extract info from rez
spikeTimes     = rez.st3(:,1);
spikeClusters  = 1+rez.st3(:,5);
spikeTemplates = rez.st3(:,2);

% get raw data around spiketimes
% WAVE = NaN(size(dat,1),numel(win),numel(spikeTimes));
WAVE = zeros(size(dat,1),numel(win),numel(spikeTimes), 'int16');
for i = 1:length(spikeTimes)
   spkwin = spikeTimes(i) + win; 
    WAVE(:,:,i) = dat(:,spkwin);
end

% find channel index with maximum amplitude template for each cluster
peakChannel = zeros(size(spikeClusters));
uClusters = unique(spikeClusters);
for c = 1:length(uClusters)
    clust     = uClusters(c);
    I         = spikeClusters == clust;
    templates =  unique(spikeTemplates(I));
    
   t = squeeze(range((rez.dWU(:,:,templates)),1));
   m = max(max(t));
   if any(size(t) ==1)
       chidx = find(t == m);
   else
       [chidx, ~] = find(t == m);
   end
   peakChannel(I) = chidx;
    
end

out.spikeTimes    = spikeTimes;
out.spikeClusters = spikeClusters;
out.spikeWaves    = WAVE;
out.peakChannel   = peakChannel;
% DEV, to add:
% template waveforms
% channels spanned by each cluster