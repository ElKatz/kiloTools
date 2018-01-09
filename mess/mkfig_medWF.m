function hF = mkfig_medWF(sp)
%   hF = mkfig_medWF(sp)
%
% makes figure with a subplot per cluster. each subplot has the median
% waveform on each of the nChannels. works on kilosort output sp struct
%
% INPUT:
%   sp  - works on kilosort output sp struct (see getSp.m)

yBuffer = 1e4;
figure,
for iClu = 1:sp.nClu
    subplot(1, sp.nClu, iClu)
    hold on
    for iCh = 1:sp.nCh
        plot(((iCh-1)*yBuffer) + squeeze(sp.medWFs(iClu, iCh,:))');
    end
    xlim([0 75])
    set(gca, 'XTick',[])
    xlabel('Time')
    ylim([-yBuffer sp.nCh*yBuffer])
    set(gca, 'YTick',[])
    if iClu==1
        ylabel('Channel')
    end
    title(['cid: ' num2str(sp.cids(iClu))])
end