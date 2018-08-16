function [hA] = plot_probeVoltage(samples)

%%
hold on;

skip = 1e2;

nCh      = size(samples,1);
nSamples = size(samples, 2);

samples = single(samples);
chStd    = std(samples(1,:));
buffer   = 15*chStd;
yTickVal = nan(nCh,1);
for iCh = 1:nCh
    plot(-iCh*buffer + samples(iCh,1:skip:end));
    yTickVal(iCh) = -iCh*buffer;
end
ylabel('channel')
xlabel('sample')
set(gca, 'YTick', flipud(yTickVal), 'YTickLabel', nCh:-1:1)
