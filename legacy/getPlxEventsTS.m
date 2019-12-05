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

