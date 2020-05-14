% ad hoc script to go over datasets that have been converted and sorted
% before I had implemented the analog data extraction that is now part of
% the conversion process.



%% first, test:
% get a plx file:
plPath = 'Y:\LAB PROJECTS\scDualVprobe\data\ram20190627a\ram20190627a.pl2';

clear ad
ai(1) = PL2Ad(plPath, 'AI01');
ai(2) = PL2Ad(plPath, 'AI02');
ai(3) = PL2Ad(plPath, 'AI03');
ai(4) = PL2Ad(plPath, 'AI04');

%
figure;
for ii = 1:4
    subplot(4,1,ii);
    plot(ai(ii).Values(1:1e5))
end

%% dir into a directory, and create an aiChannels for EACH pl file:
D = dir('Y:\LAB PROJECTS\scInactivateAndRecord\data');
idxGood = arrayfun(@(x) ~contains(x.name, '.'), D);
D = D(idxGood);
dirList = arrayfun(@(x)fullfile(x.folder, x.name), D, 'UniformOutput', 0);
disp(dirList)
nDirs = numel(dirList);
plPath = cell(nDirs,1);
trackGood = false(nDirs,1);
tic
for iD = 16:18
    % for iD = 1:1
    files       = dir(dirList{iD});
    idxPl       = arrayfun(@(x) contains(x.name, '.pl'), files);
    plPath{iD}  = fullfile(files(idxPl).folder, files(idxPl).name);
    disp('--------------------------')
    disp('getting analog data for :')
    disp(plPath{iD});
    
    
    % define output folder as 'kiloSorted':
    outFolder = fullfile(dirList{iD}, 'kiloSorted');
    
    % track whether the folder exists or not. If it doesn't then just skip
    % and ask questions later.
    if exist(outFolder, 'dir')
        trackGood(iD) = true;
    else
        trackGood(iD) = false;
        disp('oh oh. No kiloSorted folder. Weird. Skip now ask later.')
        continue;
    end
    
    % extract analog channels:
    clear ai
    ai(1) = PL2Ad(plPath{iD}, 'AI01');
    ai(2) = PL2Ad(plPath{iD}, 'AI02');
    ai(3) = PL2Ad(plPath{iD}, 'AI03');
    ai(4) = PL2Ad(plPath{iD}, 'AI04');
    
    % construct a vector of time (in seconds) that corresponds to the
    % voltages in ad.Values.
    ii = 1; % all ad channels are identical, so I'll jsut do it for 1:
    aiTimeStamps = zeros(sum(ai(ii).FragCounts),1);
    
    % sample duration
    sampDur = 1/ai(ii).ADFreq;
    
    % how many fragments of recording?
    nFrags = length(ai(ii).FragTs);
    currentSample = 1;
    for i = 1:nFrags
        chunkIndex = currentSample:(currentSample + ai(ii).FragCounts(i) - 1);
        chunkTimeStamps = ai(ii).FragTs(i) + (0:(ai(ii).FragCounts(i)-1))*sampDur;
        aiTimeStamps(chunkIndex) = chunkTimeStamps;
        currentSample = chunkIndex(end)+1;
    end
    
    %% extract info:
    pl2 = PL2GetFileIndex(plPath{iD});
    
    
    %% save
    
    save(fullfile(outFolder, 'aiChannels.mat'), 'aiTimeStamps', 'ai');
    save(fullfile(outFolder, 'info_pl2.mat'), 'pl2');
    disp(['Took me ' num2str(toc) ' seconds'])
end




%% manual versoin:
tic
dirList = 'Y:\LAB PROJECTS\scInactivateAndRecord\data\ramsey20200320
\';
files       = dir(dirList);
idxPl       = arrayfun(@(x) contains(x.name, '.pl'), files);
plPath{iD}  = fullfile(files(idxPl).folder, files(idxPl).name);
disp('--------------------------')
disp('getting analog data for :')
disp(plPath{iD});

% define output folder as 'kiloSorted':
%     outFolder = fullfile(dirList, 'kiloSorted');
    outFolder = dirList;
    % track whether the folder exists or not. If it doesn't then just skip
    % and ask questions later.
    if exist(outFolder, 'dir')
        trackGood(iD) = true;
    else
        trackGood(iD) = false;
        disp('oh oh. No kiloSorted folder. Weird. Skip now ask later.')
    end
    
    % extract analog channels:
    clear ai
    ai(1) = PL2Ad(plPath{iD}, 'AI01');
    ai(2) = PL2Ad(plPath{iD}, 'AI02');
    ai(3) = PL2Ad(plPath{iD}, 'AI03');
    ai(4) = PL2Ad(plPath{iD}, 'AI04');
    
    % construct a vector of time (in seconds) that corresponds to the
    % voltages in ad.Values.
    ii = 1; % all ad channels are identical, so I'll jsut do it for 1:
    aiTimeStamps = zeros(sum(ai(ii).FragCounts),1);
    
    % sample duration
    sampDur = 1/ai(ii).ADFreq;
    
    % how many fragments of recording?
    nFrags = length(ai(ii).FragTs);
    currentSample = 1;
    for i = 1:nFrags
        chunkIndex = currentSample:(currentSample + ai(ii).FragCounts(i) - 1);
        chunkTimeStamps = ai(ii).FragTs(i) + (0:(ai(ii).FragCounts(i)-1))*sampDur;
        aiTimeStamps(chunkIndex) = chunkTimeStamps;
        currentSample = chunkIndex(end)+1;
    end
    
    % extract info:
    pl2 = PL2GetFileIndex(plPath{iD});
    
    
    % save
    
    save(fullfile(outFolder, 'aiChannels.mat'), 'aiTimeStamps', 'ai');
    save(fullfile(outFolder, 'info_pl2.mat'), 'pl2');
    disp(['Took me ' num2str(toc) ' seconds'])
beep








