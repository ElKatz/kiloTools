function getPlxEventsTS

% choose PLX files
[f, p] = uigetfile('*.plx', 'MultiSelect', 'Off');

plxFile = [p, f];
% load file time stamp value at the start of each segment
[~,~,tsStartVals, fn, ~] = plx_ad(plxFile, 28);

% get strobed TS and vals from plx file
strobed = get_strobed(plxFile);


% We're going to construct "tsMap" that has a timestamp for every sample
% recorded. We'll do this in steps, one for each PLX file. The tsMap vector
% is used to convert from the spike index output of kiloSort (these are
% simply integers that indicate which sample number each spike occurred
% at), to a time (in seconds) relative to the beginning of the first PLX
% file recording. This is needed because the event time stamps (eventTS)
% from the PLX files are in seconds. Because the event time stamps in each
% PLX file are referenced to time in that file, we need to adjust
% successive segments of the tsMap vector and the event time stamps for the
% PLX file accordingly.

% sampling rate
fs = 4 * 10^4;

% loop over recording segments
currSample = 1;
for i = 1:length(fn)
    timeStamps = (0:(fn(i)-1))/fs;
    tsMap(currSample:(currSample - 1 + fn(i))) = timeStamps + tsStartVals(i);
    currSample = currSample + fn(i);
end

% write strobed TS and values, and tsMap to the firectory with the PLX file(s).
save([p, f(1:end-4) '.mat'], 'strobed', 'tsMap', '-v7.3');
end


function strobed = get_strobed(filename)
%% load stim codes list
codes = stimcodes_FST;
disp('Loading stimcodes_FST');
%% load Events
disp('Loading Plx event ts');
[~, ts, sv] = plx_event_ts(filename,257);
%% get strobed
disp('extracting trial event codes');
% strobed = get_strobed(sv,codes);
start = find(sv == codes.trialBegin);
stop =  find(sv == codes.trialEnd);
j = 1; strobed = {};
if (length(start) ~= length(stop)) || any((stop-start)<=0)
    disp('mismatched start and stop codes but carrying on fixing');
end
while ~isempty(start) && ~isempty(stop)
    if length(start) >=2 && length(stop) >= 1
        if start(1) < stop(1) && stop(1) < start(2)
            strobed{j} = [sv(start(1):stop(1)) ts(start(1):stop(1))];
            j = j+1;
            start(1)=[];stop(1)=[];
        elseif stop(1) < start(1) % start code missed
            stop(1) = [];
        elseif start(2) < stop(1) % stop code missed
            start(1) = [];
        end
    elseif length(start) == 1 && length(stop) == 1 && start(1) < stop(1)
        strobed{j} = [sv(start(1):stop(1)) ts(start(1):stop(1))];
        j = j+1;
        start(1)=[];stop(1)=[];
    else
        break
    end
end
end