function [hFig] = clustSummary(clust, opts)
%   [hF] = mkfig.clustSummary(clust, opts)


if ~exist('opts', 'var')
    opts = [];
end

nUnits  = numel(clust); 

% set the number of plots you intend to plot for each unit.
% this is hard coded and should match the largest 'plotNum' value below.
% e.g. 1=mean waveform, 2=isi dist, 3=spikeCount,4=text...
nPlots = 7; 

hFig            = figure;
inchesPerRow    = 2;
inchesPerCol    = 2;
figSz           = [inchesPerCol * nPlots, inchesPerRow * nUnits];
clr             = lines(nUnits);

% sort units by peak channel number:
[~, idxSort] = sort([clust.peakCh]);

% init:
iPlot = 1;

for iU = idxSort
    
    %% waveform:
    plotNum = 1;
    subplot(nUnits, nPlots, (iPlot-1)*nPlots + plotNum);
    clust(iU).plot_waveform(clr(iPlot,:));
    %% isi dist 
    plotNum = plotNum+1;
    subplot(nUnits, nPlots, (iPlot-1)*nPlots + plotNum);
    clust(iU).plot_isiDist(clr(iPlot,:));
    
    %% isi dist zoom
    plotNum = plotNum+1;
    subplot(nUnits, nPlots, (iPlot-1)*nPlots + plotNum);
    clust(iU).plot_isiDistZoom(clr(iPlot,:));
    
    %% spike count over time:
    plotNum = plotNum+1;
    subplot(nUnits, nPlots, (iPlot-1)*nPlots + plotNum);
    clust(iU).plot_spCountOverTime(clr(iPlot,:));
    
    %% ACG:
    plotNum = plotNum+1;
    subplot(nUnits, nPlots, (iPlot-1)*nPlots + plotNum);
    clust(iU).plot_acg(clr(iPlot,:));
    
     %% ACG zoom:
    plotNum = plotNum+1;
    subplot(nUnits, nPlots, (iPlot-1)*nPlots + plotNum);
    clust(iU).plot_acgZoom(clr(iPlot,:));
    
    %% text:
    plotNum = plotNum+1;
    subplot(nUnits, nPlots, (iPlot-1)*nPlots + plotNum);
    clust(iU).plot_text;
    
    iPlot = iPlot + 1;
end

supertitle(clust(1).info.dsn, 12)
formatFig(hFig, figSz, 'nature')
if isfield(opts, 'saveFigs') && opts.saveFigs == true
    if ~isfield(opts, 'dirFigs')
        opts.dirFigs = pwd;
    end
    saveas(hFig, fullfile(opts.dirFigs, 'figures', 'clustSummary.pdf'));
end
    