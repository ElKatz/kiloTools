function hF = medWfOverChannels(sp)
%   hF = mkfig.medWfOverChannels(sp)
%
% makes figure with a subplot per cluster. each subplot has the median
% waveform on each of the nChannels. works on kilosort output sp struct
%
% INPUT:
%   sp  - works on kilosort output sp struct (see getSp.m)
%%
hFig = figure;
figSz = [sp.nClusters, 4];

yScaler = .003;
x  = 1:size(sp.medWfs, 3);

for iClu = 1:sp.nClusters
    subplot(1, sp.nClusters, iClu)
    hold on
    for iCh = 1:sp.nChannels
%         plot(((iCh-1)*yBuffer) + squeeze(sp.medWFs(iClu, iCh,:))');
        plot(x + sp.xcoords(iCh), sp.ycoords(iCh) + yScaler * squeeze(sp.medWfs(iClu, iCh,:))');
    end
    xlim([-50 350])
    ylim([min(sp.ycoords) - 100, max(sp.ycoords) + 100])
    set(gca, 'XTick',[], 'YTick',[])
    if iClu==1
        ylabel('x distance (µ)')
        xlabel('y distance (µ)')
        set(gca, 'XTick', unique(sp.xcoords), 'YTick', unique(sp.ycoords))
    end
    title(['cid: ' num2str(sp.clusterId(iClu))])
end

formatFig(hFig, figSz);
saveas(hFig, 'testFig.pdf')