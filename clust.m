classdef clust < handle
    
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
        medWfOnPeakCh   = [];       % The median waveform on its peak channel
    end
    
    methods
        
        function obj = sp2clust(sp, ksDir)
            %
            % converts the session-level 'sp' struct into a structarray for single
            % untis: 'su', of size nUnits.
            
            disp('Converting ''sp'' struct into single-unit struct ''su'':...')
            % number of single units:
            nClu = numel(sp.clusterId);
           
            for iS = 1:nClu
                fprintf('unit %0.0f\r', iS)
                obj(iS).id              = sp.clusterId(iS);
                spikeIdx                = sp.spikeClusters == sp.clusterId(iS);
                obj(iS).times           = sp.spikeTimesSecs(spikeIdx);
                obj(iS).score           = sp.clusterScore(iS);
                obj(iS).uQ              = sp.uQ(iS);
                obj(iS).cR              = sp.cR(iS);
                obj(iS).isiV            = sp.isiV(iS);
                if ~isempty(sp.wf)
                    obj(iS).wf           = squeeze(sp.wf(:, :, spikeIdx));
                else
                    obj(iS).wf           = [];
                end
                obj(iS).medWfPeakCh      = sp.peakCh(iS);
                obj(iS).medWfOnPeakCh    = sp.medWfOnPeakCh(iS,:);
                %     su(iS).medWfPerChannel  = squeeze(sp.medWfs(iS,:,:))';
                %     medWfAmpPerChannel      = max(su(iS).medWfPerChannel) - min(su(iS).medWfPerChannel);
                %     su(iS).medWfPeakChannel             = find(medWfAmpPerChannel == max(medWfAmpPerChannel), 1);
                %     su(iS).medWfOnPeak     = su(iS).medWfPerChannel(:, su(iS).medWfPeakChannel);
                % %     su(iS).medWfOnPeak     = median(squeeze(su(iS).wf(peakChannel, :, :)),2);
                obj(iS).info       = sp.info;
            end
            
            disp('Done!')
            
            %% save su
            if exist('ksDir', 'var')
                save(fullfile(ksDir, 'clu.mat'), 'obj')
                disp('Done saving ''su''')
            else
                disp('no ''ksDir'' provided as input so I''m not saving nada')
            end
        end
        
    end
end