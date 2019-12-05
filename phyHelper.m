function phyHelper(IDs, varargin)

%%
p = inputParser;
p.addOptional('dir', [])
p.addOptional('task', 'tod')
p.addOptional('analysis', 'default')
p.parse(varargin{:});

dir = p.Results.dir;
task = p.Results.task;

if isempty(dir)
    dir = uigetdir();
end

%%

% where's the current master codes file:
addpath('Z:\PLDAPS_vK2_MASTER\')
codes = pds.initCodes;

% load the "stimCodes", also store some info about the stimcodes so we can
% figure out which ones we have in the file
fns         = fieldnames(codes);
codeNums    = cellfun(@(x)codes.(x),fns);
doPlots     = true;
binWidth    = 0.04;



% load spike times
sp.spikeTimesSecs = readNPY([dir filesep 'spike_times_seconds.npy']);

% load event info:
eventInfoFileName = [dir filesep 'strobedEvents.mat'];
ev = load(eventInfoFileName);
ev = ev.strobedEvents;

% load most recent cluster information
sp.spikeClusters = readNPY([dir filesep 'spike_clusters.npy']);

% in-line function definition
flat = @(x)x(:);

% build the PSP vector for convolving, it's 200ms long.
fs              = 1/1000;
pspWinDur       = 0.02;
pspFun          = @(t)(1-exp(-t/fs)).*exp(-t/pspWinDur);
pspVector       = pspFun(0:fs:pspWinDur);
pspVector       = pspVector(2:end) / sum(pspVector(2:end)); % normalize

%%
% make a cell-array consisting of all the event values ("codes") for
% each trial, starting with the "trialBegin" value.
trialBegins = [0; ev.stopTs];
nTrials     = length(ev.stopTs);
eventValues = cell(nTrials, 1);
eventTimes  = eventValues;
for i = 1:nTrials
    % find event timestamps greater than the start time of the current
    % trial and smaller than the start time of the next trial.
    tempEventTimeLogical = ev.eventInfo.Ts > trialBegins(i) & ...
        ev.eventInfo.Ts < trialBegins(i+1);
    
    % store event times and event values corresponding to the logical
    % vector defined above.
    eventValues{i}  = ev.eventInfo.Strobed(tempEventTimeLogical);
    eventTimes{i}   = ev.eventInfo.Ts(tempEventTimeLogical);
end


% task codes
taskCode = afterCode(codes, 'taskCode', eventValues);


switch task
    case 'tod'
        % find TOD
        todTrials = taskCode == codes.uniqueTaskCode_tod;
        
        % find SOA
        soaValues = afterCode(codes, 'soaDuration_x1000', eventValues);
        
        % find good TOD trials
        correctOrNot = afterCode(codes, 'correctOrnot_plus10', eventValues);
        
        % logical index to zero SOA trials
        zeroSoaTrials = todTrials & soaValues == 0 & ...
            ismember(correctOrNot, [10 11]);
        
        % time of target onset in zero SOA trials
        todTargetOnsetTime = timeInTrial(codes, 'targ1On', eventValues, ...
            eventTimes);
        
        %%%%% HERE'S WHERE YOU CHANGE THE PLOTS %%%%%%%
        gc(1, 1, :) = {zeroSoaTrials, todTargetOnsetTime, [-0.1 1.2], 'target onset', 'both targets', []};
        %     gc(1, 2, :) = {sacLoc2, targetOnsetTimes, [-0.3 0.3], '', 'target location 2', []};
        
        %     gc(2, 1, :) = {sacLoc1, saccadeOnsetTimes, [-0.1 0.5], 'saccade onset', 'target location 1', []};
        %     gc(2, 2, :) = {sacLoc2, saccadeOnsetTimes, [-0.3 0.3], '', 'target location 2', []};
    otherwise
end

nIDs = length(IDs);

for k = 1:nIDs
    
    spikeTS = sp.spikeTimesSecs(sp.spikeClusters == IDs(k));
    % parse the spikeTimes into those occurring in each trial.
    spTimes  = arrayfun(...
        @(x,y)spikeTS(spikeTS>=x & spikeTS<=y), ...
        trialBegins(1:end-1), trialBegins(2:end), ...
        'UniformOutput', false);
    
    % data holders
    if k == 1
        nr      = size(gc,1);
        nc      = size(gc,2);
        yPSP    = cell(nr, nc, nIDs);
        xPSP    = cell(nr, nc);
        meanPSP = yPSP;
        n       = yPSP;
        nws     = yPSP;
    end
    
    fh = figure('Position', [1 + (k-1)*801 400 800 400], 'Color', [1 1 1], ...
        'MenuBar', 'None', 'ToolBar', 'None');
    mc = parula(nc + 2);
    mc = mc(2:end-1,:);
    for i = 1:nr
        ax(i) = mySubPlot([1,nr,i]);
        hold on
        h = NaN(nc, 1);
        for j = 1:nc
            if ~isempty(gc{i,j,1})
                try
                    [yPSP{i,j,k}, xPSP{i,j}, meanPSP{i,j,k}, ...
                        n{i,j,k}, nws{i,j,k}] = pspOutput(...
                        gc{i,j,1}, ...
                        gc{i,j,2}, ...
                        spTimes, ...
                        pspVector, fs, gc{i,j,3}, gc{i,j,6});
                catch me
                    keyboard
                end
                
                if ~isempty(xPSP{i,j}) && ~isempty(meanPSP{i,j,k})
                    h(j) = plot(xPSP{i,j}, meanPSP{i,j,k}, ...
                        'Color', mc(j, :), 'LineWidth', 2);
                end
            end
        end
        legendObject = legend(h(ishandle(h)), gc(i, ishandle(h), 5));
        legendObject.FontSize = 12;
        legendObject.Box = 'Off';
        ylabel(gc{i,1,4}, 'FontSize', 16);
        try
            set(ax(i), 'XLim', gc{i,1,3});
        catch me
            keyboard
        end
    end
    set(gcf, 'NextPlot', 'add')
    set(ax, {'TickDir'}, {'Out'}, 'NextPlot', 'add');
    arrayfun(@(x)plot(x, [0 0], ...
        [0 max(flat(get(x, 'YLim')))], 'k--'), ax);
    if ~isempty(which('suptitle'))
        suptitle(['cluster ID = ' num2str(IDs(k))]);
    else
        supertitle(['cluster ID = ' num2str(IDs(k))]);
    end
    
%     close(fh);
end
end

function [y,x] = pspIfy(spikeTimes, pspVector, fs, startWin, endWin, tt)

% start and end times (either the first and last spike, or the beginning /
% end of the window depending on which is "larger").
startTime   = fix(min([spikeTimes(:); startWin])/fs)*fs;
endTime     = fix(max([spikeTimes(:); endWin])/fs)*fs;

% this vector will be zero everywhere except samples in which a spike has
% occured.
X           = zeros(round((endTime - startTime)/fs)+1, 1);

if isempty(spikeTimes)
    y = NaN(length(X), 1);
else
    % set sample-times in X when spikes occured (relative to "startTime") to 1
    si          = fix((spikeTimes-startTime)/fs)+1;
    X(si)       = 1;
    
    % convolve spiketrain with psp-vector, only keep the portion that's
    % the same length as "X".
    y           = conv(X, pspVector, 'same')*(1/fs);
    
    % make a time index for the part we actually care abev.
    timeIndex   = startTime:fs:endTime;
    
    % make a logical index for the part we care about
    g           = timeIndex-startWin > -fs/2 & timeIndex-endWin < fs/2;
    
    % It's possible the first/last spike in the trial occurs in the window of
    % interest. In this case, we want to replace everything before the first
    % and after the last spike with NaNs. But NOT if there's NO spikes in the
    % window. To do this, construct a logical index of all the time points
    % that THIS CAUSES WEIRD THINGS TO HAPPEN. STOP DOING IT.
    % gNaN        = (timeIndex - min(spikeTimes) > -fs/2) & (timeIndex - max(spikeTimes) < fs/2);
    
    % if nnz(gNaN) ~= length(y)
    %   y(~gNaN)     = NaN;
    % end
    
    % set y beyond "tt" to NaNs
    y(timeIndex > tt) = NaN;
    
    % limit the convolution output to the length of X
    y           = y(1:length(X));
    
    % limit the convolution output further to the window of interest
    y           = y(g);
    
    % limit the time vector to the window of interest
    x           = timeIndex(g);
    
    if abs(min(x)-startWin)>(fs/2) || abs(max(x)-endWin)>(fs/2)
        keyboard
    end
end

end

function out    = afterCode(codes, field, eventValues)

% loop over the cell array "eventValues" and find the eventValue that
% follows the first occurrance of "codes.(field)" in each array element.
nEventValues = length(eventValues);
for i = 1:nEventValues
    tempFieldIndex = find(eventValues{i} == codes.(field));
    if ~isempty(tempFieldIndex)
        if length(eventValues{i}) >= tempFieldIndex + 1
            temp = eventValues{i}(tempFieldIndex + 1);
            out(i,1) = temp(1);
        else
            out(i,1) = 0;
        end
    else
        out(i,1) = 0;
    end
end

end

function y                          = funINE(fun, args)

if iscell(args)
    emptyBool = any(cellfun(@isempty, args));
else
    emptyBool = isempty(args);
end
if ~emptyBool
    y = fun(args{:});
else
    y = NaN;
end

end

function time   = timeInTrial(codes, field, eventValues, eventTimes)
time            = cellfun(@(x,y)x(y==codes.(field)), eventTimes, eventValues, 'UniformOutput', false);
end

function [yPSP, xPSP, meanPSP, n, nws]      = pspOutput(group, time, spikeTimes, pspVector, fs, xlimin, varargin)

% if a "truncTime" argument was passed, use that to exclude activity for
% spikes occuring after "truncTime". Otherwise, set "truncTime" to a
% "large" value so that spikes aren't unintentionally excluded.
if ~isempty(varargin{1})
    truncTime = varargin{1};
else
    truncTime = num2cell(NaN(length(time),1));
    truncTime(group) = cellfun(@(x)x + 100, time(group), ...
        'UniformOutput', false);
end

% calculate "truncTimeAligned": relative to whatever event we're aligning
% spikes to.
truncTimeAligned = cellfun(@(x,y)x-y, truncTime(group), time(group), ...
    'UniformOutput', false);

% calculate "timeAlignedSpikes": spike times relative to whatever event
% we're aligning to.
timeAlignedSpikes = cellfun(@(x,y)x-y, spikeTimes(group), time(group), ...
    'UniformOutput', false);

% psp-ify all the spike-trains, excluding spikes that are greater than
% "truncTime".
yPSP = cellfun( ...
    @(x, y)pspIfy(x, pspVector, fs, xlimin(1), xlimin(2), y), ...                         % anonymous function definition
    timeAlignedSpikes, ... % time-aligned spikes
    truncTimeAligned, ... % exclusion time
    'UniformOutput', false);

xPSP = xlimin(1):fs:xlimin(2);
meanPSP = nanmean(cell2mat(yPSP'),2);
n   = nnz(group);
nws = nnz(~cellfun(@isempty, spikeTimes(group)));

end