function [spikeTimes] = thresholdSpikes(time, voltages, stdMultiplier, nChannels)

spikeTimes = cell(nChannels, 1);
tic;
for iCh = 1:nChannels
    v = double(voltages(iCh,:));
    distLimits = prctile(v(1:1e3:end), [.05 99.95]); % take 99% of the (subsampled) voltage distribution
    idx = (v < distLimits(1) | v > distLimits(2));
    v(idx) = nan;
    fprintf('nan''d %d, values\n', sum(idx))
    sd = nanstd(v);
    thresh = (-stdMultiplier * sd);
    assert(sign(thresh)==-1, 'your threshold is positive. that''s a no-no');
    belowThreshIdx = v < thresh;
    spikeIdx = [0 diff(belowThreshIdx)] == 1;
    spikeTimes{iCh} = time(spikeIdx);
end
toc

%% 
% 
% figure, hold on
% plot(time(1:1e2:end), v(1:1e2:end));
% line(xlim, [thresh thresh], 'Color', 'r')
% 
