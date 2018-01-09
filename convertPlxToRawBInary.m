function [samples, datPath] = convertPlxToRawBInary(fullPathPlx, opts)
%   [samples, datPath] = convertPlxToRawBInary(fullPathPlx, opts)
%
% Converts continuous data from plx file to matrix 'samples' of size:
% [nSamples, nChannels] and saves as binary .dat file in same folder as plx
% file (unless specified otherwise in opts)
% INPUT:
%   fullPathPlx - full path to plx file that contains continuous data, or
%                 path to folder that contains mulitple plx files to be
%                 concatenated.
%                   if left empty, function opens
%
%   opts - (optional) struct of options:
%   .outputFolder - full path to folder to save dat file (default: same
%                   folder as plx file)
%
% OUTPUT:
%   samples - [nChannels, nSamples] consisting of all continuous data
%   datPath - path to the .dat file

% 20170417 - lnk & arb: bare bones
% 20170525 - lnk: added concatenation option
%                 added saving of strobes, time stamps
%                 added mat bkp

% 2do:
% need to figure out how to get continuous channel string. is it the same
% for amar's plx files and my mpx coverted to plx files?


%% paths:

addpath(genpath('~/Dropbox/Code/spike_sorting/toolboxes/Matlab Offline Files SDK'));
addpath(genpath('~/Dropbox/Code/repos/pdstools/dependencies'));

%% get plx files:
if ~exist('fullPathPlx', 'var')
    [plxFileName, plxFolder] = uigetfile('*.plx', 'Select plx file/s for conversion', '~/Dropbox/Code/spike_sorting/pilot_datasets/', 'MultiSelect', 'on');
    %     fullPathPlx = '/Users/leorkatz/Dropbox/Code/spike_sorting/pilot_datasets/amar_vProbe_dirTun_170510/dir_170510_2.plx';
    % fullPathPlx = 'C:\EPHYS on SSD\katz\test\amar_recording_probe\170413_mrg-01.plx';
    % fullPathPlx = 'C:\EPHYS on SSD\amar\rawdata\170413_1_mrg-01.plx';
    % fullPathPlx = 'C:\EPHYS on SSD\katz\F170320-0001Post.plx';
else
    [plxFolder, plxFileName, ext]  = fileparts(fullPathPlx);
    plxFileName = [plxFileName ext];
end
disp('--------------------------------------------------------------')
disp('converting plx file/s from:')
disp(plxFolder)
disp('--------------------------------------------------------------')

%% optional arguements:
% init:
if ~exist('opts', 'var')
    opts = struct;
end

% output folder:
if exist('opts', 'var') &&  isfield(opts, 'outputFolder')
    outputFolder = opts.outputFolder;
else
    outputFolder = plxFolder;
end

%% get file names:
if ~iscell(plxFileName)
    plxFileName = mat2cell(plxFileName, 1);
end
disp(plxFileName')

% raw binary .dat file name: named after the first plx file:
datPath = fullfile(outputFolder, [plxFileName{1} '.dat']);

%% begin conversion:

nFiles = numel(plxFileName);
for iF = 1:nFiles
    disp(['performing conversion on:   ' plxFileName{iF}])
    
    % full path to plx file:
    fullPathPlx = fullfile(plxFolder, plxFileName{iF});
    
    % create a list of all ad continuous channel names in cell array:
    [nCh, adChName]     = plx_adchan_names(fullPathPlx);
    chNameList = cell(nCh,1);
    for ii = 1:nCh
        chNameList{ii} = adChName(ii,:);
    end
    
    % get indices for the spike channels (i.e. not the lfp):
    % *alphaLab spike files are saved with the name "SPK"  and plexon are
    % "CSPK" so I'm just gonna use "SPK" as the identifier:
    spkChannelStr = 'SPK';
    idxSpkCh      = false(numel(chNameList),1);
    for iCh = 1:numel(chNameList)
        if strfind(chNameList{iCh}, spkChannelStr)
            idxSpkCh(iCh) = true;
        else
            idxSpkCh(iCh) = false;
        end
    end
    
    % get number of spikes counts per ad channel and get indices for those that
    % have data:
    [~, samplecounts] = plx_adchan_samplecounts(fullPathPlx);
    idxDataCh = samplecounts~=0;
    
    % get indices for channels that are both spk channels & have data:
    idxGoodCh =  idxSpkCh & idxDataCh;
    
    % nChannels & nSamples:
    nChannels   = sum(idxGoodCh);
    tmp         = samplecounts(idxGoodCh);
    nSamples    = tmp(1); % taking the number of samples in first spk channel. Rest are identical.
    
    % build data matrix 'samples' of size [nChannels, nSamples]:
    samples     = zeros(nChannels, nSamples, 'int16');
    tChRead     = nan(nChannels,1); % time keeping
    
    % gotta map out indices to plxeon's ad channel numbers:
    [~,   adChNumber]   = plx_ad_chanmap(fullPathPlx);
    spkChNumber = adChNumber(idxGoodCh);
    
    disp(['Getting data from ' num2str(sum(idxGoodCh)) ' spike channels!'])
    for iCh = 1:nChannels
        tic
        [adfreq, n, ts, fn, ad] = plx_ad(fullPathPlx, spkChNumber(iCh)); % retunrs signal in miliVolts
        if n>0
            tChRead(iCh) = toc;
            fprintf('Took me %0.3f secs to read channel #%0.0d \n', tChRead(iCh), spkChNumber(iCh));
            % data matrix 'samples':
            samples(iCh,1:n) = int16(ad);
        else
            disp('wtf?!')
            keyboard
        end
    end
    
    
    
    % %
    % %     % identify the continuous channels by name:
    % %     [nCh, adChName]     = plx_adchan_names(fullPathPlx);
    % %     [~,   adChNumber]   = plx_ad_chanmap(fullPathPlx);
    % %     chNameList = cell(nCh,1);
    % %     for ii = 1:nCh
    % %         chNameList{ii} = adChName(ii,:);
    % %     end
    % %     % continuous spike channels are identified by the strong below. This is
    % %     % not a great method but it works for alphaLab files (which are named
    % %     % "SPK") and plexon recorded files (named "CSPK")
    % %     continuousChannelStr = 'SPK';
    % %     % get indices for continuous channels:
    % %     idx         = false(numel(chNameList),1);
    % %     for iCh = 1:numel(chNameList)
    % %         if strfind(chNameList{iCh}, continuousChannelStr)
    % %             idx(iCh) = true;
    % %         else
    % %             idx(iCh) = false;
    % %         end
    % %     end
    % %     contChNumber = adChNumber(idx); % adding 1 because adChNumber is zero-based list
    % %
    % %     % make sure that the channel numbers are similar between files:
    % %     if iF==1
    % %         contChNumberFirstFile = contChNumber;
    % %     else
    % %         assert(mean(contChNumber==contChNumberFirstFile)==1, 'goddamn! There''s a mismatch between the continuous channel numbers of this file vs 1st file. wtf?');
    % %     end
    % %
    % %     % build data matrix 'samples' of size [nChannels, nSamples]:
    % %     nChannels   = numel(contChNumber);
    % %     [~, ~, ~, contCounts] = plx_info(fullPathPlx, 1);
    % %     nSamples    = contCounts(contChNumber(1));
    % %     samples     = zeros(nChannels, nSamples, 'int16');
    % %     tChRead     = nan(nChannels,1);
    % %     chGood      = false(nChannels,1);
    % %     hWait = waitbar(0, 'Loading up channel by channel and populating ''samples''');
    % %     for iCh = 1:nChannels
    % %         tic
    % %         [adfreq, n, ts, fn, ad] = plx_ad(fullPathPlx, contChNumber(iCh)); % retunrs signal in miliVolts
    % %         if n>0
    % %             chGood(iCh) = true;
    % %             tChRead(iCh) = toc;
    % %             sprintf('Took me %0.3f secs to read channel %0.0d', tChRead(iCh), contChNumber(iCh));
    % %
    % %             % data matrix 'samples':
    % %             samples(iCh,1:n) = int16(ad);
    % %         else
    % %             chGood(iCh) = false;
    % %         end
    % %         waitbar(iCh/nChannels, hWait)
    % %     end
    % %     close(hWait)
    % %     disp(['Got data from' num2str(sum(chGood)) ' continuous channels!'])
    % %
    % %     % take only good channels:
    % %     samples = samples(chGood, :);
    
    % append data to .dat file:
    fidout = fopen(datPath, 'a'); % opening file for appending
    fwrite(fidout, samples, 'int16');
    fclose(fidout);
    
    %     fSamples{iF} = samples; for debugging
end % iF


%% TEST ZONE
% clear t
% for iCh = 1:nChannels
%     tic;
%     pl = readPLXFileC(fullPathPlx, 'continuous', continuousChannelNumbers(iCh));
%     t.singleChannelLoad_readPlx(iCh) = toc;
%
%     tic;
%     [~,~,~,~, ad] = plx_ad(fullPathPlx, continuousChannelNumbers(iCh));
%     t.singleChannelLoad_plx_ad(iCh) = toc;
% end

%% look at activity
nSecs = 3;
% fs = 40000;
% figure, hold on
% for iCh = 1:size(samples,2)
%     plot(iCh*500 + samples(1:10:nSecs*fs, iCh))
% end
% set(gca, 'XTick', 0:fs/2:nSecs, 'XTickLabel', 0:.5:nSecs)

%
%
%
% toc