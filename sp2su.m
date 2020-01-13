function su = sp2su(sp, ksDir)
%   su = sp2su(sp, ksDir)
%
% converts the session-level 'sp' struct into a structarray for single
% untis: 'su', of size nUnits.
% This is mostly for convenience as 'su' holds unit-specific information
% (eg spike times, ISI violations, unit score etc...)
%
% INPUT:
%   sp - the 'sp' struct created by getSp.m
% 
% OUTPUT:
%   su - the 'su' stuct with unit-specific data. 


disp('Converting ''sp'' struct into single-unit struct ''su'':...')
% number of single units:
nSus = numel(sp.clusterId);

su  = struct;

for iS = 1:nSus
    fprintf('unit %0.0f\r', iS)
    su(iS).clusterId        = sp.clusterId(iS);
    spikeIdx                = sp.spikeClusters == sp.clusterId(iS);
    su(iS).times            = sp.spikeTimesSecs(spikeIdx);
    su(iS).clusterScore     = sp.clusterScore(iS);
    su(iS).uQ               = sp.uQ(iS);
    su(iS).cR               = sp.cR(iS);
    su(iS).sp.isiV_fpRate   = sp.sp.isiV_fpRate(iS);
    su(iS).sp.isiV_rate     = sp.sp.isiV_rate(iS);
    if ~isempty(sp.wf)
        su(iS).wf           = squeeze(sp.wf(:, :, spikeIdx));
    else
        su(iS).wf           = [];
    end
    su(iS).peakCh           = sp.peakCh(iS);
    su(iS).medWfOnPeakCh    = sp.medWfOnPeakCh(iS,:);
%     su(iS).medWfPerChannel  = squeeze(sp.medWfs(iS,:,:))';
%     medWfAmpPerChannel      = max(su(iS).medWfPerChannel) - min(su(iS).medWfPerChannel);
%     su(iS).medWfPeakChannel             = find(medWfAmpPerChannel == max(medWfAmpPerChannel), 1);
%     su(iS).medWfOnPeak     = su(iS).medWfPerChannel(:, su(iS).medWfPeakChannel);
% %     su(iS).medWfOnPeak     = median(squeeze(su(iS).wf(peakChannel, :, :)),2);
    [a,b,c]             = fileparts(sp.info.ksDir);
    su(iS).info.dsn     = b;
    su(iS).info.ksDir   = ksDir; 
    su(iS).info.Fs      = sp.sample_rate;
    su(iS).info.meta    = 'kiloSort';
end

disp('Done!')

%% save su
if exist('ksDir', 'var')
    save(fullfile(ksDir, 'su.mat'), 'su')
    disp('Done saving ''su''')
else
    disp('no ''ksDir'' provided as input so I''m not saving nada')
end