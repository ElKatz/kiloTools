% test:
% get a plx file:
filename = 'Y:\LAB PROJECTS\scDualVprobe\data\ram20190627a\ram20190627a.pl2';

clear ad
ad(1) = PL2Ad(filename, 'AI01');
ad(2) = PL2Ad(filename, 'AI02');
ad(3) = PL2Ad(filename, 'AI03');
ad(4) = PL2Ad(filename, 'AI04');

%%
figure;
for ii = 1:4
    subplot(4,1,ii);
    plot(ad(ii).Values(1:1e5))
end

%%
D = dir('Y:\LAB PROJECTS\scDualVprobe\data');
idxGood = arrayfun(@(x) ~contains(x.name, '.'), D);
D = D(idxGood);
dirList = arrayfun(@(x)fullfile(x.folder, x.name), D, 'UniformOutput', 0);
nDirs = numel(dirList);
plPath = cell(nDirs,1);
trackGood = false(nDirs,1);
tic
for iD = 1:nDirs
% for iD = 1:1
    files       = dir(dirList{iD});
    idxPl       = arrayfun(@(x) contains(x.name, 'pl'), files);
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
    
    % extract info:
    pl2 = PL2GetFileIndex(plPath{iD});
    
    % extract analog channels:
    clear ad
    ad(1) = PL2Ad(plPath{iD}, 'AI01');
    ad(2) = PL2Ad(plPath{iD}, 'AI02');
%     ad(3) = PL2Ad(plPath{iD}, 'AI03');
%     ad(4) = PL2Ad(plPath{iD}, 'AI04');
    
    save(fullfile(outFolder, 'plAnalogInChannels'), 'pl2', 'ad');
    disp(['Took me ' num2str(toc) ' seconds'])
end