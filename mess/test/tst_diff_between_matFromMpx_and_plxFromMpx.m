% compare the matFromMpx to plxFromMpx

%% matFromMpx 

rawPath = '/Users/leorkatz/Dropbox/Code/fst_attention/data/matFromMpx/snap/20170502/raw.mat';
a = load(rawPath);

TS{1} = a.tEvent;
SV{1} = a.vEvent;

spk_TS{1} = spikeTimes{1};
N_SPK{1} = numel(spk_TS{1});

figure, 
subplot(221); hold on
plot(TS{1})
plot(spk_TS{1})

subplot(222); 
hist(SV{1}, 50)



%% plxFromMpx
datapath_mat = '/Users/leorkatz/Dropbox/Code/fst_attention/data/plx_and_mat/snap/20170502/F170502-0004_mrg.mat';
q = load(datapath_mat);
TS{2} = q.eventTime_mrg;
SV{2} = q.eventValue_mrg;

datapath_plx = '/Users/leorkatz/Dropbox/Code/fst_attention/data/plx_and_mat/snap/20170502/F170502-0004 - CopyPost_mrg-01.plx';
pl  = readPLXFileC(datapath_plx, 'spikes');
[N_SPK{2}, spk_TS{2}] = plx_ts(datapath_plx, pl.SpikeChannels(idx).Channel, 0);


subplot(223); hold on
plot(TS{2})
plot(spk_TS{2})

subplot(224); 
hist(SV{2}, 50)

%%
ii=1;

ts = TS{ii};
sv = SV{ii};
spk_ts = spk_TS{ii};
n_spk = N_SPK{ii};
