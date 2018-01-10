function mpx = enforceEqualSmplingRateAcrossChannels(mpx, chString)
%
% In alphaLab, it may be the case that some spk channels were recorded with
% some sampling rate and others with another. This function sets the lowers
% common sampling rate for all channels:

%%

% only perform this for spk channels, not lfp:
if isempty(strfind(lower(chString), 'spk'))
    return;
end

%% get channel names:

chNameList  = {mpx.Channel_ID_Name_Map.Name};

% get indices for the spike channels (i.e. not the lfp):
idxSpkCh = getIdx_forChannelType(chNameList, 'SPK');

idxGoodCh = idxSpkCh; % should I check for anything else? is there an 'idxGood' to make? check if all are zero otr someting? apply theo zero xhexk to plx too?

% good spk channles names:
spkChName   = chNameList(idxGoodCh);

%% Equate sampling rate:


nChannels = sum(idxGoodCh);
% nSamples: should be identical accross channels but lets check:
nSamplesCh = nan(nChannels,1);
for iCh = 1:nChannels
    nSamplesCh(iCh) = numel(mpx.(spkChName{iCh}).Samples);
end
if numel(unique(nSamplesCh))==1
    % all good, we have a single sampling rate
    fprintf('\n\tAll spk channels recorded at similar sampling rate');
    fprintf('\n\tNo need to enforce equal sampling');
    return;
else
    fprintf('\n\tEnforcing equal sampling rate across spk channels...')
    nSampsLo = min(nSamplesCh);
    ptrFsLo = find(nSamplesCh == nSampsLo,1);
    fsLo = mpx.(spkChName{ptrFsLo}).KHz;
    ptrFsHi = find(nSamplesCh ~= nSampsLo);
    for iHi = ptrFsHi
        fprintf('\n\tboo. ch #%0.0d was recorded at %0.0d Khz while the min Fs is %0.0d. Subsampling...', iHi, mpx.(spkChName{iHi}).KHz, fsLo)
        nSampsHi = numel(mpx.(spkChName{iHi}).Samples);
        subSampleInteger = nSampsHi ./ nSampsLo;
        
        if abs(subSampleInteger - round(subSampleInteger)) < 1e-3
            % all good. there a 1 sample difference in one of the channels.
            subSampleInteger = round(subSampleInteger);
        else
            error('must be an integer otherwise I''m fucked');
        end

        mpx.(spkChName{iHi}).Samples = mpx.(spkChName{iHi}).Samples(1:subSampleInteger:end);
        mpx.(spkChName{iHi}).KHz = mpx.(spkChName{iHi}).KHz / subSampleInteger;
    end
end


