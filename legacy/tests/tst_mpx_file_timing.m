mpx = load('/Users/iMacKatz/Dropbox/_transfer_big_data/20170617/f170617-0005.mat');

%%
spk     = mpx.CSPK_001;
tSpk    = linspace(mpx.CSPK_001_TimeBegin, mpx.CSPK_001_TimeEnd, numel(spk));

% event smaple, time, & value (s, t, & v):
sEvent  = mpx.CInPort_001(1,:);
tEvent  = mpx.CInPort_001(1,:) ./ (mpx.CInPort_001_KHz * 1e3);
vEvent  = mpx.CInPort_001(2,:);

%% plot spikes actvity and event time:

figure(11), clf; hold on
plot(tSpk(1:1e3:end), spk(1:1e3:end))

% all events:
% plot(tEvent, 1e-3.*bsxfun(@times, vEvent, ones(1,numel(vEvent))), '.');

% just trial start:
ptr         = find(vEvent==30001); %trial start
visScaler   = 15;
plot(tEvent(ptr), visScaler .* ones(1, numel(ptr)), '.', 'Color', [.2 .7 .2], 'MarkerSize', 6)
% just trial end:
ptr         = find(vEvent==30009); %trial end
visScaler   = 15;
plot(tEvent(ptr), visScaler .* ones(1, numel(ptr)), '.', 'Color', [1 .2 .2], 'MarkerSize', 3)
%% plot ttl pulses:

% find the sample at which TTL channel begins:
sTtlStart   = find(tSpk >= mpx.CTTL_001_TimeBegin, 1);
sUp         = mpx.CTTL_001_Up + sTtlStart;
sDn         = mpx.CTTL_001_Down + sTtlStart;

% build logical vector of good/bad samples:
idxGood     = false(numel(tSpk),1); % init

try
    % mpx files are inconsistent.  sometimes there are more sUp, other 
    % times there are more sDn. gotta cater to each file's individual needs..
    sDn(1) = [];
    sDn(end+1) = numel(tSpk);
    for ii = 1:numel(sUp)
        idxGood(sUp(ii):sDn(ii)) = true;
    end
catch me
    error('this didn''t work due to idiodyncracies in this particualr mpx file. Look at the number of sUp vs. sDn and fix pairing.')
end
visScaler = 14;
plot(tSpk, visScaler .* idxGood)
%%
%
%
%% check number of sUp vs. sDn in a file:
fileList = dir;
for iF = 1:numel(fileList)
    fullPath = fullfile(fileList(iF).folder, fileList(iF).name);
    [pth, name, ext] = fileparts(fullPath);
    if strcmp(ext, '.mat')
        load(fullPath)
        disp(fullPath)
        [size(CTTL_001_Up)' size(CTTL_001_Down)']
    end
end
