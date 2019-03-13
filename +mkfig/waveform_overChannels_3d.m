function hFig = waveform_overChannels_3d(sp, opts)
%   hFig = mkfig.waveform_overChannels_3d(sp, opts)

if ~exist('opts', 'var')
    opts = [];
end

%%
hFig = figure; 
figSz = [8 8];

plot_waveform_overChannels_3d(sp);

%%
formatFig(hFig, figSz);
supertitle(sp.info.dsn, 12)
formatFig(hF, figSz)
if isfield(opts, 'saveFigs') && opts.saveFigs == true
    if ~isfield(opts, 'dirFigs')
        opts.dirFigs = pwd;
    end
    saveas(hF, fullfile(opts.dirFigs, 'figures', 'unitSummary.pdf'));
end