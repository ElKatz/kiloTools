function hA = plot_waveform_overChannels(sp)
hold on
title('Waveform')

chList          = 1:sp.n_channels_dat;
chListFlip      = fliplr(chList);
chListFlipNeg   = -chListFlip;
yScaler         = 0.003;

for iS = 1:sp.nClusters
    ch = sp.medWfPeakCh(iS);
    wf = sp.medWfOnPeakCh(iS,:) .* yScaler;
    wfLength = numel(wf);
    x  = (1:wfLength)  - ceil(wfLength/2);
    hP(iS) = plot(x, -ch + wf); % flipping sign for plotting
    xJitter = 0;
    yJitter = rand()*2-1;
    text(-wfLength/2 + xJitter, -ch + yJitter, ['cid ' num2str(sp.clusterId(iS))], 'Color', hP(iS).Color, 'FontSize', 6)
end 
ylim([chListFlipNeg(1)-1 chListFlipNeg(end)+1])
xlim(.75.*[-wfLength wfLength])
ylabel('Channel Number')
set(gca, 'YTick', chListFlipNeg, 'YTickLabel', chListFlip)
