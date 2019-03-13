function hA = plot_spikeCount_overChannles(sp)
% hA = plot_spikeCount_overChannles(sp)
%%

title('Spike Count')
chList          = 1:sp.n_channels_dat;
chListFlip      = fliplr(chList);
chListNeg       = -chList;
chListFlipNeg   = fliplr(chListNeg);
spCountPerCh    = nan(sp.n_channels_dat,10); % designing the array this way allows for spCounts from multiple units to be later plotted as 'stacked'

for iS = 1:sp.nClusters
    ch = sp.medWfPeakCh(iS);
    spCount         = sum(sp.spikeTimesSamps(sp.spikeClusters == sp.clusterId(iS)));
    if spCount > 0
        % replace 1st nan with spCount:
        ptr = find(isnan(spCountPerCh(ch,:)),1);
        spCountPerCh(ch,ptr) = spCount;
    end
end 
hB = barh(chListNeg, spCountPerCh, 'stacked');
ylabel('Channel Number')
set(gca, 'YTick', chListFlipNeg, 'YTickLabel', chListFlip)

