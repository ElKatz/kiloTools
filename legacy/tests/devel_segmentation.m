%% test

% plotting the strobe/event times vs the ttl times

% q = load('/Users/iMacKatz/Dropbox/Code/spike_sorting/pilot_datasets/leor_vProbe_attn_20170511_cont/f170511-0019.mat');
mpx = load('/Users/iMacKatz/Dropbox/Code/spike_sorting/pilot_datasets/leor_vProbe_attn_20170511_cont/f170511-0004.mat');

visBool = false;
%% Get the spk channels

% first makes sure that all channels are recorded at the same smapling
% rate. If not, then go with the lowest common denominator.

% get indices for the spike channels (i.e. not the lfp):
idxSpkCh = get_idxForSpkChannels(chName);

idxGoodCh = idxSpkCh; % should I check for anything else? ios there an 'idxGood' to make? check if all are zero otr someting? apply theo zero xhexk to plx too?

% take note of the good spk channles for posterity:
spkChName   = chName(idxGoodCh);

% nChannels:
nChannels   = numel(spkChName);
% nSamples: should be identical accross channels but lets check:
nSamplesCh = nan(nChannels,1);
for iCh = 1:nChannels
    nSamplesCh(iCh) = numel(mpx.(spkChName{iCh}).Samples);
end
if numel(unique(nSamplesCh))==1
    % all good, we have a single sampling rate
else
    nSampsLo = min(nSamplesCh);
    ptrFsLo = find(nSamplesCh == nSampsLo,1);
    fsLo = mpx.(spkChName{ptrFsLo}).KHz;
    ptrFsHi = find(nSamplesCh ~= nSampsLo);
    for iHi = ptrFsHi
        disp(['boo. ch #' num2str(iHi) ' was recorded at ' num2str(mpx.(spkChName{iHi}).KHz) ' Khz while the min Fs is ' num2str(fsLo) '. Subsampling...'])
        subSampleInteger = numel(mpx.(spkChName{iHi}).Samples) ./ nSampsLo;
        assert(mod(subSampleInteger,1)==0, 'error: must be an integer otherwise I''m fucked');
        mpx.(spkChName{iHi}).Samples = mpx.(spkChName{iHi}).Samples(1:subSampleInteger:end);
        mpx.(spkChName{iHi}).KHz = mpx.(spkChName{iHi}).KHz / subSampleInteger;
    end
end
% Now that they are all the saem Fs, taking the number of
% samples in first spk channel.
nSamples = numel(mpx.(spkChName{1}).Samples);

% build data matrix 'samples' of size [nChannels, nSamples]:
samples     = zeros(nChannels, nSamples, 'int16');
%             tChRead     = nan(nChannels,1); % time keeping

disp(['Getting data from ' num2str(nChannels) ' spike channels!'])
for iCh = 1:nChannels
    tic
    samples(iCh,:) = mpx.(spkChName{iCh}).Samples;
end

%% get event times:
% event samples, time & values (s, t, v):
% sEv = q.CInPort_001.Samples(1,:);
tEv = mpx.CInPort_001.Samples(1,:) ./ (mpx.CInPort_001.KHz * 1e3);
vEv = mpx.CInPort_001.Samples(2,:);

% external trigger up & down, sample and time::
% sUp = q.CTTL_001.Up;
% sDn = q.CTTL_001.Down;
tUp = mpx.CTTL_001.Up ./ (mpx.CTTL_001.KHz * 1e3) + mpx.CTTL_001.TimeBegin;
tDn = mpx.CTTL_001.Down ./ (mpx.CTTL_001.KHz * 1e3)+ mpx.CTTL_001.TimeBegin;

% could use the strobed start/stop trial instead of TTL for consistency


if visBool
    % makse sure that times match between events and ttl up/dn
    figure,
    subplot(211); hold on; title('tEv')
    hist(tEv,50)
    subplot(212); hold on; title('tUp')
    hist(tUp,50)
end

%%

% for every channel, create time vector:
tSpk    = linspace(mpx.CSPK_001.TimeBegin, mpx.CSPK_001.TimeEnd, numel(mpx.CSPK_001.Samples));
idxGood = false(numel(tSpk,1)); % init

% if the file started with the ttl already up there will not be a tUp that
% corresponds to the tDn. Giving tUp a 0 to account for that:
if tDn(1) < tUp(1)
    tUp = [0 tUp];
end
% if the file ended with the ttl already up then there will not be a
% corresponding Dn. Giving tDn the largest element in tSpk
if tUp(end) > tDn(end)
    tDn = [tDn tSpk(end)];
end
assert(numel(tUp)==numel(tDn), 'error: different number of elements in tUp and tdn')
assert(all(sign(tDn - tUp))==1, 'error: not all tDn come after correspoinding tUp')

nTrigger = numel(tUp);
ptrUp = nan(nTrigger, 1);
ptrDn = nan(nTrigger, 1);
for ii = 1:nTrigger
    ptrUp(ii) = find(tSpk >= tUp(ii),1);
    ptrDn(ii) = find(tSpk >= tDn(ii),1);
    idxGood(ptrUp(ii):ptrDn(ii)) = true;
end

if visBool
    % visualize
    % makse sure that up/dn signals correspond to the spike times
    % makse sure that events are at trial times and not in between
    figure, hold on
    plot(tSpk, mpx.CSPK_001.Samples)
    plot(tSpk, (idxGood-.5)*100)
    plot(tEv, 60*ones(numel(tEv,1)), 'o')
end

%% SEGMENT:

samples     = samples(:, idxGood);
tSamples    = tSpk(idxGood);

if visBool
    figure, plot(tSamples, samples(1,:))
end
