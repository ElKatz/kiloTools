function hFig = waveformOverChannels_perCluster(sp, opts)
%   hFig = mkfig.waveformOverChannels_perCluster(sp, opts)
%
% makes figure with a subplot per cluster. each subplot has the median
% waveform on each of the nChannels. works on kilosort output sp struct
%
% INPUT:
%   sp  - works on kilosort output sp struct (see getSp.m)
%%

if ~exist('opts', 'var')
    opts = [];
end

%%
hFig = figure;
figSz = [sp.nClusters, 4];

yScaler = .1;
x  = 1:size(sp.medWfs, 3);

nX = numel(unique(sp.xcoords));
nY = numel(sp.ycoords);

[~, idxSort] = sort(sp.medWfPeakCh);

iPlot = 1;
for iClu = idxSort
    subplot(1, sp.nClusters, iPlot)
    hold on
    for iCh = 1:nY
%         plot(((iCh-1)*yBuffer) + squeeze(sp.medWFs(iClu, iCh,:))');
        plot(x + sp.xcoords(iCh) - 0.5*length(x), sp.ycoords(iCh) + yScaler * squeeze(sp.medWfs(iClu, iCh,:))');
    end
    xlim([min(sp.xcoords)-100 max(sp.xcoords)+100])
    ylim([min(sp.ycoords) - 100, max(sp.ycoords) + 100])
    set(gca, 'XTick',[], 'YTick',[])
    if iPlot==1
        ylabel('y distance (µ)')
        xlabel('x distance (µ)')
        set(gca, 'XTick', unique(sp.xcoords), 'YTick', unique(sp.ycoords))
    end
    title(['cid: ' num2str(sp.clusterId(iClu))])
    iPlot = iPlot + 1;
end

%%
formatFig(hFig, figSz);
supertitle(sp.info.dsn, 12)
if isfield(opts, 'saveFigs') && opts.saveFigs == true
    if ~isfield(opts, 'dirFigs')
        opts.dirFigs = pwd;
    end
    saveas(hFig, fullfile(opts.dirFigs, 'figures', 'waveformOverChannels_perCluster.pdf'));
end

