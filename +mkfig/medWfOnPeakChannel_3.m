function hFig = medWfOnPeakChannel_3(sp)
%   hFig = mkfig.medWfOnPeakChannel_3(sp)

%%
hFig = figure; 
figSz = [8 8];
%%
hold on

depth = -25000;
dTip2lowestCh = 500 + 800; % tip to lowest + lowest to ch28

yScaler = 0.001;

for iS = 1:sp.nClusters
    ch = sp.medWfPeakCh(iS);
    wf = sp.medWfOnPeakCh(iS,:) .* yScaler;
    x  = 1:numel(wf);
    
    % 1st z
    z   = 1.*ones(1,numel(wf)); 
    plot3(sp.xcoords(ch) + x, z, (depth + dTip2lowestCh) + sp.ycoords(ch) + wf)
    
    % 2nd z
    z   = 1.05.*ones(1,numel(wf)); 
    plot3(sp.xcoords(ch) + x, z, (depth + dTip2lowestCh) + sp.ycoords(ch) + wf)
    
    % 2nd z
    z   = 1.1.*ones(1,numel(wf)); 
    plot3(sp.xcoords(ch) + x, z, (depth + dTip2lowestCh) + sp.ycoords(ch) + wf)
end 
view(3)
% ylim([depth -15000])
% xlim([-500 500])
xlabel('x distance (µ)')
ylabel('y distance from bottom of grid (µ)')
formatFig(hFig, figSz);