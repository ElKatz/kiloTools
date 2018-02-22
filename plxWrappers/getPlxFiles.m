function plxFiles = getPlxFiles(plxDir, plxString)

% AI
% FP
% RSTART
% RSTOP
% SPK
% SPKC
% Strobed

if ~exist('plxString', 'var')
    plxString = '.pl2';
end

fileList    = dir(plxDir);
fileNames   = {fileList.name};
idx         = cellfun(@(x) strfind(x, str), fileNames, 'UniformOutput', 0);
plxFiles    = fileNames(idx);

