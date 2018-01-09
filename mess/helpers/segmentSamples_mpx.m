function [segSamples, segSamplesTime] = segmentSamples_mpx(samples, mpx, chString)
%
% function takes the 'samples' matrix [nChannels x nSamples] and segments
% it according to that start/stop TTL pulses that are sent at the beginning
% (TTL goes up) and end (TTL goes down) of every trial.

%%
% visulization bolean:
visBool = false;


%% create time vector from ch1:
if strfind(lower(chString), 'spk')
    chName = 'CSPK_001';
elseif strfind(lower(chString), 'lfp')
    chName = 'CLFP_001';
end
tSamples    = linspace(mpx.(chName).TimeBegin, mpx.(chName).TimeEnd, numel(mpx.(chName).Samples));
idxGood     = false(numel(tSamples,1)); % init

%% external trigger: 

% external trigger up & down, sample ('s') and time ('t'):
sUp = mpx.CTTL_001.Up;
sDn = mpx.CTTL_001.Down;
if sUp(1)==0 
    sUp(1) = 1; % ie first sample. 0 makes no sense
end
tUp = sUp ./ (mpx.CTTL_001.KHz * 1e3) + mpx.CTTL_001.TimeBegin;
tDn = sDn ./ (mpx.CTTL_001.KHz * 1e3) + mpx.CTTL_001.TimeBegin;   

% if the file started with the ttl already up there will not be a tUp that
% corresponds to the tDn. Giving tUp a 0 to account for that:
if tDn(1) < tUp(1)
    tUp = [0 tUp];
    sUp = [1 sUp];
end
% if the file ended with the ttl already up then there will not be a
% corresponding Dn. Giving tDn the largest element in tSpk
if tUp(end) > tDn(end)
    tDn = [tDn tSamples(end)];
    sDn = [sDn numel(tSamples)];
end
assert(numel(tUp)==numel(tDn), 'error: different number of elements in tUp and tdn')
assert(all(sign(tDn - tUp))==1, 'error: not all tDn come after correspoinding tUp')

% could use the strobed start/stop trial instead of TTL for consistency


if vis;ool
    % make sure that times match between events and ttl up/dn
    % get event times:

    tEv = mpx.CInPort_001.Samples(1,:) ./ (mpx.CInPort_001.KHz * 1e3);
    vEv = mpx.CInPort_001.Samples(2,:);

    figure,
    subplot(211); hold on; title('tEv')
    hist(tEv,50)
    subplot(212); hold on; title('tUp')
    hist(tUp,50)
end

%%


nTrigger = numel(tUp);
ptrUp = nan(nTrigger, 1);
ptrDn = nan(nTrigger, 1);
for ii = 1:nTrigger
    ptrUp(ii) = find(tSamples >= tUp(ii),1);
    ptrDn(ii) = find(tSamples >= tDn(ii),1);
    idxGood(ptrUp(ii):ptrDn(ii)) = true;
end

if visBool
    % visualize
    % makse sure that up/dn signals correspond to the spike times
    % makse sure that events are at trial times and not in between
    figure, hold on
    plot(tSamples(1:1e3:end), mpx.(chName).Samples(1:1e3:end))
    plot(tSamples, (idxGood-.5)*100)
    plot(tEv, 60*ones(numel(tEv,1)), 'o')
end

%% SEGMENT:

segSamples     = samples(:, idxGood);
segSamplesTime = tSamples(idxGood);

if visBool
    figure, plot(segSamplesTime, segSamples(1,:))
end
