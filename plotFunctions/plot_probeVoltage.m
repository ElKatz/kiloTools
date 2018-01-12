function [hA] = plot_probeVoltage(samples)

%%
hold on;

nCh      = size(samples,1);
nSamples = size(samples, 2);

chStd    = std(samples(1,:));
buffer   = 10*chStd;
yTickVal = nan(nCh,1);
for iCh = 1:nCh
    plot(-iCh*buffer + samples(iCh,:));
    yTickVal(iCh) = -iCh*buffer;
end
ylabel('channel')
xlabel('sample')
set(gca, 'YTick', flipud(yTickVal), 'YTickLabel', nCh:-1:1)
