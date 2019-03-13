function hFig = waveform_overChannels(sp, opts)
%   hFig = mkfig.waveform_overChannels(sp, opts)

if ~exist('opts', 'var')
    opts = [];
end

%%
hFig = figure; 
figSz = [3 8];

plot_waveform_overChannels(sp);

%%
supertitle(sp.info.dsn, 12)
formatFig(hFig, figSz)
if isfield(opts, 'saveFigs') && opts.saveFigs == true
    if ~isfield(opts, 'dirFigs')
        opts.dirFigs = pwd;
    end
    saveas(hFig, fullfile(opts.dirFigs, 'figures', 'waveform_overChannels.pdf'));
end
