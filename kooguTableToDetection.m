function t = kooguTableToDetection(ravenFile,soundFolder,siteCode,classification)
%t = ravenTableToStruct(ravenFile,soundFolder,siteCode,classification)
% Load a Raven Selection table of annotations, and convert this into a 
% Matlab struct, 'a', with some basic metadata for an acoustic detection.
% The metadata in Raven Selection Tables is highly variable, and here we
% assume that each Selection Table contains the following columns/fields:
%   Begin Date Time, 
%   Delta Time (s), 
%   Low Freq (Hz), 
%   High Freq (Hz)
% To make the detection data structure a bit more useful, we optionally add
% a siteCode to identify the recording, optionally link each detection to a
% soundFolder to facilitate loading acoustic data, and optionally add a
% classification code.
if nargin < 4 
    classification = '';
end
if nargin < 3
    siteCode = '';
end
if nargin < 2
    soundFolder = '';
end
warning('off','MATLAB:table:ModifiedVarnames');
if istable(ravenFile)
    t = ravenFile;
else
    t = readtable(ravenFile,'delimiter','\t');
end

nDetect = height(t);
if (nDetect == 0)
    return;
end
%% Create XBAT-like fields to facilitate compatibility with XBAT 'events'
if strcmp(t.Properties.VariableNames,'Channel')
    t.channel = t.Channel;
else
    t.channel = ones(height(t),1);
end

try
    t0 = zeros(height(t),1);
    % Assume wavFolderInfo has been called correctly and cache exists
    wavInfo = wavFolderInfo(soundFolder);
    parfor i = 1:height(t)
        ix = contains({wavInfo.fname},t.BeginFile(i));
        t0(i) = wavInfo(ix).startDate+t.FileOffset_s_(i)/86400;
    end
    t.t0 = t0;
catch
    keyboard
end
t.DeltaTime_s_ = t.EndTime_s_ - t.BeginTime_s_;
t.tEnd = t.t0+t.DeltaTime_s_/86400;        % Matlab datenum
t.duration = (t.tEnd-t.t0)*86400;          % Duration in seconds

% Older versions of Koogu
if any(strcmpi(t.Properties.VariableNames,'LowFrequency_Hz_'))
    t.fLow = t.LowFrequency_Hz_;
else % Newer versions of Koogu
    t.fLow = t.LowFreq_Hz_;
end

if any(strcmpi(t.Properties.VariableNames,'HighFrequency_Hz_'))
    t.fHigh= t.HighFrequency_Hz_;
else
    t.fHigh= t.HighFreq_Hz_;
end
t.freq = [t.fLow t.fHigh];  % Frequency vector in Hz
t.soundFolder =  cellstr(repmat(soundFolder,nDetect,1));             % Location of audio files for this detection
t.siteCode = cellstr(repmat(siteCode,nDetect,1));                    % Append site to data structure
t.classification = cellstr(repmat(classification,nDetect,1));        % Append classification to data structure;
t = sortrows(t,'t0');

