function hFig = medWfOnPeakChannel(sp)
%   hFig = mkfig.medWfOnPeakChannel(sp)

%%
hFig = figure; 
figSz = [8 8];

hold on

depth = -25000;
dTip2lowestCh = 500 + 800; % tip to lowest + lowest to ch28

yScaler = 0.03;

for iS = 1:sp.nClusters
    ch = sp.medWfPeakCh(iS);
    wf = sp.medWfOnPeakCh(iS,:) .* yScaler;
    x  = 1:numel(wf);

    plot(sp.xcoords(ch) + x, (depth + dTip2lowestCh) + sp.ycoords(ch) + wf)
end 
ylim([depth -15000])
xlim([-500 500])
xlabel('x distance (µ)')
ylabel('y distance from bottom of grid (µ)')
formatFig(hFig, figSz);