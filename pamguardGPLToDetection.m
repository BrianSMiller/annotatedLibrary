function [t, gpl] = pamguardGPLToDetection(binaryFolder, ...
    fileMask, soundFolder, siteCode, classification)
% Depends on pgmatlab from:
%   https://github.com/PAMGuard/PAMGuardMatlab/releases/tag/V1.0.0
addpath('c:\analysis\PAMGuard\pgmatlab\')
if nargin < 1 % Test dataset
    binaryFolder= 's:/data/longTermRecordings/detections/pamguard-GPS-Bp20/kerguelen2015/'
    fileMask='GPL_Detector_Generalised_Power_Law_Detector_GPL_Detections_*.pgdf';
    soundFolder = 'L:/annotatedLibrary/SORP/wav/kerguelen2015/';
    siteCode = 'kerguelen2015';
    classification = 'Bp-20';
end
gpl = loadPamguardBinaryFolder(binaryFolder,fileMask);
wavInfo = wavFolderInfo(soundFolder);
wavStartDateBins = [wavInfo.startDate wavInfo(end).startDate+1/24];
[~,filename,ext] = fileparts({wavInfo.fname}')
filename = strcat(filename,ext)
nDetect = length(gpl);
t0 = [gpl.date]';
t = table(t0);
t.channel = [gpl.channelMap]';
t.duration = [gpl.millisDuration]'/1000;% Duration in seconds
t.DeltaTime_s_ = [gpl.millisDuration]'/1000;% Duration in seconds
t.tEnd = t.t0+t.duration/86400;         % Matlab datenum

% Find file index;
[~,~,ix] = histcounts([gpl.date],wavStartDateBins)

t.BeginFile = filename(ix);
t.FileOffset_s_  = (t.t0-[wavInfo(ix).startDate]')/86400;


t.freq = cell2mat({gpl.freqLimits;}');  % Frequency vector in Hz
t.fLow = t.freq(:,1);
t.fHigh = t.freq(:,2);
t.soundFolder =  cellstr(repmat(soundFolder,nDetect,1));             % Location of audio files for this detection
t.siteCode = cellstr(repmat(siteCode,nDetect,1));                    % Append site to data structure
t.classification = cellstr(repmat(classification,nDetect,1));        % Append classification to data structure;
t = sortrows(t,'t0');