function [hFig] = unitSummary(su, opts)
%   [hF] = mkfig.unitSummary(su, opts)

if ~exist('opts', 'var')
    opts = [];
end

%% User setup:
% do you wish to plot individual waveforms behind the mean waveform?
plotWavesBool = true;
% set up the number of individual waveforms to plot behind the mean waveform
nWavesToPlot = 50;
% set number of bins within which to compute spike count:
nBinsForSpikeCount  = 100;
% set the max isi you're interested in plotting:
isiDistMax          = 0.05;
% set the number of plots you intend to plot for each unit.
% this is hard coded and should match the largest 'plotNum' value below.
% e.g. 1=mean waveform, 2=isi dist, 3=spikeCount,4=text...
nPlots = 6; 

%% figure setup:
% every unit gets a row in the figure:
nUnits  = numel(su); 

hFig            = figure;
inchesPerRow    = 2;
inchesPerCol    = 2;
figSz           = [inchesPerCol * nPlots, inchesPerRow * nUnits];
clr             = lines(nUnits);


[~, idxSort] = sort([su.medWfPeakCh]);

%%
iPlot = 1;
for iU = idxSort
    thisUnit = su(iU);
    
    %% waveform:
    plotNum = 1;
    subplot(nUnits, nPlots, (iPlot-1)*nPlots + plotNum);
    if iPlot==1
        title('med WF')
    end
    hold on
    if plotWavesBool
        if ~isempty(thisUnit.wf)
            plot(squeeze(thisUnit.wf(thisUnit.peakChannel,:, round(linspace(1,numel(thisUnit.times), nWavesToPlot)))), 'LineWidth', .5, 'Color', 'k');
        else
            disp('no individual waveforms, so skipping their plotting')
        end
    end
    
    if isfield(thisUnit, 'medWfOnPeakCh')
        plot(thisUnit.medWfOnPeakCh, 'LineWidth', 2, 'Color', clr(iU,:))
    end
%     ylim([-5e4 5e4])
    ylabel(['Cell ' num2str(thisUnit.clusterId)], 'Color', clr(iU,:))
    
    %% isi dist 
    plotNum = plotNum+1;
    subplot(nUnits, nPlots, (iPlot-1)*nPlots + plotNum);
    if iPlot==1
        title('ISI dist')
    end
    hold on
    isi = diff(thisUnit.times);
    isi(isi > isiDistMax) = [];
    histogram(isi, 50, 'FaceColor', clr(iU,:));
    xlim([0 isiDistMax])
    
    %% spike count over time:
    
    plotNum = plotNum+1;
    % this is hacky. I am not accounting for lapses in the recording. This
    % measure is only good enough for comparing units recorded at the same 
    % time, but not across different recoridngs.
    subplot(nUnits, nPlots, (iPlot-1)*nPlots + plotNum);
    if iPlot==1
        title('spCount')
    end
    hold on
    edges       = linspace(thisUnit.times(1), thisUnit.times(end), nBinsForSpikeCount+1);
    spikeCounts = histc(thisUnit.times, edges);
    spikeCounts = spikeCounts(1:end-1);
    plot(spikeCounts, 'Color', clr(iU,:), 'LineWidth', 2)
    xlim([0, nBinsForSpikeCount])
    
    %% ACG:
    plotNum = plotNum+1;
    subplot(nUnits, nPlots, (iPlot-1)*nPlots + plotNum);
    if iPlot==1
        title('ACG')
    end
    hold on
    binSize = 0.02;
    acgBins = 1/thisUnit.info.Fs/2:binSize:0.5;
    [n,x] = histdiff(thisUnit.times, thisUnit.times, acgBins);
    n = n./binSize./numel(thisUnit.times);
    stairs([-x(end:-1:1) x]-binSize/2,[n(end:-1:1) n], 'Color', clr(iU,:), 'LineWidth', 2.0);
    yl = ylim(); 
    ylim([0 yl(2)]);
    
    %% ACG zoom:
    plotNum = plotNum+1;
    subplot(nUnits, nPlots, (iPlot-1)*nPlots + plotNum);
    if iPlot==1
        title('ACG zoom')
    end
    hold on
    binRange = 0.02;
    binSize = 0.00025;
    acgBins = 1/thisUnit.info.Fs/2:binSize:binRange;
    [n,x] = histdiff(thisUnit.times, thisUnit.times, acgBins);
    n = n./binSize./numel(thisUnit.times);
    stairs([-x(end:-1:1) x]-binSize/2,[n(end:-1:1) n], 'Color', clr(iU,:), 'LineWidth', 2.0);
    hold on;
    yl2 = ylim();
    ylim([0 max([yl(2) yl2(2)])]);
    yl2 = ylim();
    plot([-0.0015 -0.0015], yl2, 'k--');
    plot([0.0015 0.0015], yl2, 'k--');
    
    %% text:
    plotNum = plotNum+1;
    subplot(nUnits, nPlots, (iPlot-1)*nPlots + plotNum);
    if iPlot==1
        title('info')
    end
    hold on
    axis([0 1 0 1])
    yPosStart = 0.1;
    if isfield(thisUnit, 'times')
        text(0.1, yPosStart, ['nSpikes = ', sprintf('%.0f', numel(thisUnit.times))]); 
        yPosStart = yPosStart+0.2;
    end
    if isfield(thisUnit, 'uQ')
        text(0.1, yPosStart, ['uQ = ', sprintf('%.2f', thisUnit.uQ)]); 
        yPosStart = yPosStart+0.2;
    end
    if isfield(thisUnit, 'cR')
        text(0.1, 0.7, ['cR = ', sprintf('%.3f', thisUnit.cR)]);
        yPosStart = yPosStart+0.2;
    end
    if isfield(thisUnit, 'isiV')
        text(0.1, 0.5, ['isiV = ', sprintf('%.3f', thisUnit.isiV)]);
        yPosStart = yPosStart+0.2;
    end
    if isfield(thisUnit, 'isolation')
        text(0.1, 0.5, ['isolation = ', sprintf('%.2f', thisUnit.isolation)]);
        yPosStart = yPosStart+0.2;
    end
    
    iPlot = iPlot + 1;
end

supertitle(su(1).info.dsn, 12)
formatFig(hFig, figSz)
if isfield(opts, 'saveFigs') && opts.saveFigs == true
    if ~isfield(opts, 'dirFigs')
        opts.dirFigs = pwd;
    end
    saveas(hFig, fullfile(opts.dirFigs, 'figures', 'unitSummary.pdf'));
end
    
    


