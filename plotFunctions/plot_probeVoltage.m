function [hA] = plot_probeVoltage(samples, skip)
%   [hA] = plot_probeVoltage(samples, skip)
%
% Quick function to plot multi-channel probe voltage traces. 
% INPUT:
%   samples - [nCh, nSamples] 'samples' matrix used in convertToDat.m 
%   skip    - optional. determines how many samples to skip in plotting
% OUTPUT:
%   hA      - handle to figure axis

%%

disp('Plotting probe voltages...')
hA = gca;
hold on;

if ~exist('skip', 'var')
    skip = 1e2;
end

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
