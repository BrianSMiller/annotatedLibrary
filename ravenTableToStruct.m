function a = ravenTableToStruct(ravenFile,soundFolder,siteCode,classification)
%a = ravenTableToStruct(ravenFile,soundFolder,siteCode,classification)
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
    ravenTable = ravenFile;
else
    ravenTable = table2struct(readtable(ravenFile,'delimiter','\t'),'toScalar',true);
end
nDetect = length(ravenTable.Channel);     % Raven selection tables will always have a channel column

%% Create XBAT-like fields to facilitate compatibility with XBAT 'events'
a.channel = ravenTable.Channel;                             % Raven channels start at 1, just like Matlab
a.t0 = datenum(ravenTable.BeginDateTime,...
                       'yyyy/mm/dd HH:MM:SS.fff');                  % Matlab datenum
a.tEnd = a.t0+ravenTable.DeltaTime_s_/86400;        % Matlab datenum
a.duration = (a.tEnd-a.t0)*86400;           % Duration in seconds
a.fLow = ravenTable.LowFreq_Hz_;
a.fHigh= ravenTable.HighFreq_Hz_;
a.freq = [ravenTable.LowFreq_Hz_ ravenTable.HighFreq_Hz_];  % Frequency vector in Hz
a.soundFolder =  repmat(soundFolder,nDetect,1);             % Location of audio files for this detection
a.siteCode = repmat(siteCode,nDetect,1);                    % Append site to data structure
a.classification = repmat(classification,nDetect,1);        % Append classification to data structure;
a.fileName = ravenTable.BeginFile;
a = orderfields(a);
a = structofarrays2arrayofstructs(a);