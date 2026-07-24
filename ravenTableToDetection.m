function ravenTable = ravenTableToDetection(ravenFile,soundFolder,siteCode,classification)
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
    classification = '';   % sentinel: per-row Tags passthrough below if present, else 'none'
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
    ravenTable = readtable(ravenFile,'delimiter','\t');
end

% Remove 'From column' since this is inconsistently present
if any(strcmpi(ravenTable.Properties.VariableNames,'From'))
    ravenTable.('From') = [];
end

nDetect = height(ravenTable);
if (nDetect == 0)
    return;
end
%% Create XBAT-like fields to facilitate compatibility with XBAT 'events'
if strcmp(ravenTable.Properties.VariableNames,'Channel')
    ravenTable.channel = ravenTable.Channel;
else
    ravenTable.channel = ones(height(ravenTable),1);
end

try
ravenTable.t0 = datenum(ravenTable.BeginDateTime,...
                       'yyyy/mm/dd HH:MM:SS.fff');   % Matlab datenum
catch
    sampleRate = mean(ravenTable.BegFileSamp_samples_./ravenTable.FileOffset_s_);
    for i = 1:height(ravenTable)
        ravenTable.t0(i) = filenameToTimeStamp(ravenTable.BeginFile{i},...
            '_yyyy-mm-dd_HH-MM-SS.') +...
            ravenTable.BegFileSamp_samples_(i)/sampleRate /86400;
    end
end
if ~any(strcmpi(ravenTable.Properties.VariableNames,'DeltaTime_s_'))
    ravenTable.DeltaTime_s_ = ravenTable.EndTime_s_-ravenTable.BeginTime_s_;
end
ravenTable.tEnd = ravenTable.t0+ravenTable.DeltaTime_s_/86400;        % Matlab datenum
ravenTable.duration = (ravenTable.tEnd-ravenTable.t0)*86400;          % Duration in seconds
ravenTable.fLow = ravenTable.LowFreq_Hz_;
ravenTable.fHigh= ravenTable.HighFreq_Hz_;
ravenTable.freq = [ravenTable.LowFreq_Hz_ ravenTable.HighFreq_Hz_];  % Frequency vector in Hz
ravenTable.soundFolder =  cellstr(repmat(soundFolder,nDetect,1));             % Location of audio files for this detection
ravenTable.siteCode = cellstr(repmat(siteCode,nDetect,1));                    % Append site to data structure

% classification: explicit value (unchanged behaviour) if given, otherwise
% per-row Tags passthrough when a Tags column exists, otherwise 'none'.
% FIXED 2026-07-24: repmat(classification,nDetect,1) with the default ''
% stayed 0x0 regardless of nDetect (repmat of an empty array doesn't grow),
% so cellstr(...) produced a 0-row result -- a row-count mismatch on
% assignment for ANY nDetect>0 whenever this function was called without
% an explicit 4th argument (e.g. via matchbox-style 3-arg converter calls).
if isempty(classification)
    if any(strcmp('Tags', ravenTable.Properties.VariableNames))
        ravenTable.classification = cellstr(string(ravenTable.Tags));      % per-row passthrough
    else
        ravenTable.classification = repmat({'none'}, nDetect, 1);          % no Tags column to fall back on
    end
else
    ravenTable.classification = cellstr(repmat(classification,nDetect,1)); % explicit value, unchanged behaviour
end

end
