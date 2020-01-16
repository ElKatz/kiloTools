classdef clustClass < handle
    % An individual cluster-based class that allows for single cluster
    % based analyses.
    %
    % clust = clustClass;
    %   Initializes the clust class.
    %
    % clust.sp2clust_forClusterId(sp, clusterId);
    %   Converts the 'sp' struct into the 'clustClass' with all the
    %   properties listed below.
    %   'sp' - output of function getSp.m (from kiloTools github repo)
    %   'clusterId' -
    %
    
    properties(GetAccess='public', SetAccess='protected')
        info            = [];       % Can't hurt to carry some info
        id              = [];       % cluster ID
        times           = [];       % spike times (s)
        score           = [];       % 0=noise, 1=mua, 2=good
        uQ              = [];       % unit Quality (see sqKilosort.computeAllMeasures.m)
        cR              = [];       % contamination ratio (see sqKilosort.computeAllMeasures.m)
        isiV_rate       = [];       % isi violation rate (see sqKilosort.computeAllMeasures.m)
        isiV_fpRate     = [];       % isi violation false positive rate (see sqKilosort.computeAllMeasures.m)
        wf              = [];       % individual waveforms
        peakCh          = [];       % the channel which has the peak amplitude waveform
        distToTip       = [];       % distnace from peakCh to electrode tip (in order to compute depth given tip depth)
        medWfOnPeakCh   = [];       % The median waveform on its peak channel
        defaultClr      = [];       % default color for this unit in plots
    end
    
    methods(Static=true)
        
        function obj = clustClass(sp, clusterId)
            
            if nargin < 2
                error('Must provide ''sp'' struct and ''clusterId'' for which you wish to constrcut this clustClass');
            end
            disp(['Constructing clustClass for clusterId ' num2str(clusterId)]);
            
            % this cluster ID is unit# 'iS':
            iS = find(sp.clusterId == clusterId);
            obj.id              = sp.clusterId(iS);
            spikeIdx            = sp.spikeClusters == sp.clusterId(iS);
            obj.times           = sp.spikeTimesSecs(spikeIdx);
            obj.score           = sp.clusterScore(iS);
            obj.uQ              = sp.uQ(iS);
            obj.cR              = sp.cR(iS);
            if isfield(sp, 'isiV_fpRate')
                obj.isiV_fpRate     = sp.isiV_fpRate(iS);
                obj.isiV_rate       = sp.isiV_rate(iS);
            else
                [obj.isiV_fpRate, numViolations] = ISIViolations(obj.times, 0.0005, 0.001);
                nSpikes         = length(obj.times);    
                obj.isiV_rate   = numViolations/nSpikes;
            end
            if ~isempty(sp.wf)
                obj.wf          = squeeze(sp.wf(:, :, spikeIdx));
            else
                obj.wf          = [];
            end
            if isfield(sp, 'peakCh')
                obj.peakCh          = sp.peakCh(iS);
            else
                obj.peakCh          = sp.medWfPeakCh(iS);
            end
            obj.medWfOnPeakCh   = sp.medWfOnPeakCh(iS,:);
            obj.info            = sp.info;
            obj.info.Fs         = sp.sample_rate;
            tmpClr              = lines(64);
            obj.defaultClr      = tmpClr(iS);
        end
    end
    methods 
        %%
        function plot_waveform(obj, clr)
            nWavesToPlot = 50;
            
            title('med WF')
            hold on
            if ~isempty(obj.wf)
                plot(squeeze(obj.wf(obj.peakCh,:, round(linspace(1,numel(obj.times), nWavesToPlot)))), 'LineWidth', .5, 'Color', 'k');
            end
            xvals = (1:length(obj.medWfOnPeakCh)) / obj.info.Fs;
            plot(xvals, obj.medWfOnPeakCh, 'LineWidth', 2, 'Color', clr)
            xlim([xvals(1) xvals(end)])
            ylabel('uV')
            xlabel('Time (s)')
        end
        
        %%
        function plot_isiDist(obj, clr)
            isiDistMax          = 0.1;
            title('isi')
            hold on
            isi = diff(obj.times);
            isi(isi > isiDistMax) = [];
            histogram(isi, 30, 'FaceColor', clr);
            xlim([0 isiDistMax])
            line([0.001 0.001], ylim, 'Color', 'k', 'LineStyle', '--')
            xlabel('Time (s)')
            ylabel('Count')
        end
        
        %%
        function plot_isiDistZoom(obj, clr)
            isiDistMax          = 0.01;
            title('isi ZOOM')
            hold on
            isi = diff(obj.times);
            isi(isi > isiDistMax) = [];
            histogram(isi, 30, 'FaceColor', clr);
            xlim([0 isiDistMax])
            line([0.001 0.001], ylim, 'Color', 'k', 'LineStyle', '--')
            xlabel('Time (s)')
            ylabel('Count')
        end
        
        %%
        function plot_spCountOverTime(obj, clr)
            
            %% spike count over time:
            nBinsForSpikeCount  = 20;
            
            % this is hacky. I am not accounting for lapses in the recording. This
            % measure is only good enough for comparing units recorded at the same
            % time, but not across different recoridngs.
            
            title('Stability')
            hold on
            edges       = linspace(obj.times(1), obj.times(end), nBinsForSpikeCount+1);
            spikeCounts = histc(obj.times, edges);
            spikeCounts = spikeCounts(1:end-1);
            stairs(spikeCounts, 'Color', clr, 'LineWidth', 2)
            xlim([0, nBinsForSpikeCount])
            xlabel('bin #')
            ylabel('spCount per bin')
            
        end
        
        function plot_acg(obj, clr)
            if ~exist('clr', 'var')
                clr = lines(1);
            end
            %% ACG:
            title('ACG')
            hold on
            binSize = 0.01;
            acgBins = 1/obj.info.Fs/2:binSize:0.5;
            [n,x] = histdiff(obj.times, obj.times, acgBins);
            n = n./binSize./numel(obj.times);
            stairs([-x(end:-1:1) x]-binSize/2,[n(end:-1:1) n], 'Color', clr, 'LineWidth', 2.0);
            xlabel('Time (s)')
            ylabel('Count')
        end
        
        function plot_acgZoom(obj, clr)
            %% ACG zoom:
            title('ACG zoom')
            hold on
            binRange = 0.01;
            binSize = 0.00025;
            acgBins = 1/obj.info.Fs/2:binSize:binRange;
            [n,x] = histdiff(obj.times, obj.times, acgBins);
            n = n./binSize./numel(obj.times);
            stairs([-x(end:-1:1) x]-binSize/2,[n(end:-1:1) n], 'Color', clr, 'LineWidth', 2.0);
            hold on;
            plot([-0.001 -0.001], ylim, 'k--');
            plot([0.001 0.001], ylim, 'k--');
            xlabel('Time (s)')
        end
        
        function plot_text(obj, textSize)
            %% text:
            
            % ADD TO TEXT WHETHER CLUSTER WAS CONDIFERED GOOD OR MUA
            title('info')
            hold on
            set(gca, 'Visible', 'off')
            axis([0 1 0 1])
            if ~exist('textSize', 'var')
                textSize = 6;
            end
            text(0.05, .9, ['ID = ', sprintf('%.0f', obj.id)], 'fontSize', textSize);
            text(0.05, .8, ['score = ', sprintf('%.0f', obj.score)], 'fontSize', textSize);
            text(0.05, .7, ['nSpikes = ', sprintf('%.0f', numel(obj.times))], 'fontSize', textSize);
            text(0.05, .6, ['uQ = ', sprintf('%.2f', obj.uQ)], 'fontSize', textSize);
            text(0.05, .5, ['cR = ', sprintf('%.3f', obj.cR)], 'fontSize', textSize);
            text(0.05, .4, ['isiV rate = ', sprintf('%.3f', obj.isiV_rate)], 'fontSize', textSize);
            
        end
    end
end

%%
