function [idxSpkCh] = getIdx_forChannelType(chNameList, channelString, verbose)
% get indices for the spike channels (i.e. not the lfp):
% - alphaLab spike files are saved with the name "CSPK"
% - plexon spike files are saved with the name "SPK"
% so I'm just gonna use "SPK" as the identifier:

if ~exist('channelString', 'var')
    channelString = 'SPK';
end
if ~exist('verbose', 'var')
    verbose = false;
end




idxSpkCh      = false(numel(chNameList),1);
for iCh = 1:numel(chNameList)
    if strfind(chNameList{iCh}, channelString)
        idxSpkCh(iCh) = true;
    else
        idxSpkCh(iCh) = false;
    end
end

if verbose
    if sum(idxSpkCh)==0
        warning('Found 0 channels that amtch your channelString');
    else
        fprintf('Found %d channels that match your criteria\n', sum(idxSpkCh))
    end
end