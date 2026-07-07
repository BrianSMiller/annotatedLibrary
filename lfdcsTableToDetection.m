function t = lfdcsTableToDetection(ravenFile,soundFolder,siteCode,classification)
% Convert LFDCS detections into a detection table compatible with the
% callDensity workflow.
%
% Optional inputs:
%   soundFolder
%   siteCode
%   classification

%% Optional inputs
if nargin < 4
    classification = '';
end
if nargin < 3
    siteCode = '';
end
if nargin < 2
    soundFolder = '';
end

%% Read table
warning('off','MATLAB:table:ModifiedVarnames');
if istable(ravenFile)
    t = ravenFile;
else
    t = readtable(ravenFile,'delimiter','\t');
end

nDetect = height(t);
if nDetect == 0
    return;
end

%% Channel
if any(strcmpi(t.Properties.VariableNames,'channel'))
    t.channel = t.channel;
else
    t.channel = ones(nDetect,1);
end

%% Convert datetime columns
% Two supported source formats:
%   1) start_datetime already present as 'yyyy-MM-dd HH:mm:ss.SSS'
%   2) Detection_DateTime_Raw ('MM/dd/yy HH:mm:ss') + Fractional_Second,
%      as produced by the raw LFDCS export (after header-stripping via
%      standardizeLFDCSOutput.m). Reassembled here into start_datetime.
if any(strcmpi(t.Properties.VariableNames,'start_datetime'))
    t.start_datetime = datetime(t.start_datetime, ...
        'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
elseif any(strcmpi(t.Properties.VariableNames,'Detection_DateTime_Raw'))
    dt = datetime(t.Detection_DateTime_Raw, ...
        'InputFormat', 'MM/dd/yy HH:mm:ss');
    if any(strcmpi(t.Properties.VariableNames,'Fractional_Second'))
        dt = dt + seconds(t.Fractional_Second);
    end
    dt.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
    t.start_datetime = dt;
else
    error(['Could not find a start datetime column (expected ', ...
        'start_datetime or Detection_DateTime_Raw + Fractional_Second).']);
end

%% Resolve Dynamic Column Names (Flexible Mapping)
% Check for duration column variations
if any(strcmpi(t.Properties.VariableNames, 'detection_duration_s'))
    t.detection_duration_s = t.detection_duration_s;
elseif any(strcmpi(t.Properties.VariableNames, 'Duration_s'))
    t.detection_duration_s = t.Duration_s;
else
    error('Could not find a duration column (expected detection_duration_s or Duration_s).');
end

% Check for minimum frequency variations
if any(strcmpi(t.Properties.VariableNames, 'minimum_frequency_hz'))
    t.minimum_frequency_hz = t.minimum_frequency_hz;
elseif any(strcmpi(t.Properties.VariableNames, 'Min_Frequency_Hz'))
    t.minimum_frequency_hz = t.Min_Frequency_Hz;
else
    error('Could not find a minimum frequency column.');
end

% Check for maximum frequency variations
if any(strcmpi(t.Properties.VariableNames, 'maximum_frequency_hz'))
    t.maximum_frequency_hz = t.maximum_frequency_hz;
elseif any(strcmpi(t.Properties.VariableNames, 'Max_Frequency_Hz'))
    t.maximum_frequency_hz = t.Max_Frequency_Hz;
else
    error('Could not find a maximum frequency column.');
end

%% Calculate end_datetime using the dynamically mapped duration
t.end_datetime = t.start_datetime + seconds(t.detection_duration_s);

%% Handle missing file_start_time to avoid relative calculation crashes
if ~any(strcmpi(t.Properties.VariableNames,'file_start_time'))
    t.file_start_time = t.start_datetime;
else
    t.file_start_time = datetime(t.file_start_time, ...
        'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
end

%% Relative times within wav file
t.BeginTime_s_ = seconds(t.start_datetime - t.file_start_time);
t.EndTime_s_ = t.BeginTime_s_ + t.detection_duration_s;
t.DeltaTime_s_ = t.detection_duration_s;

%% Matlab datenums
t.t0 = datenum(t.start_datetime);
t.tEnd = datenum(t.end_datetime);
t.duration = t.detection_duration_s;

%% Frequency fields
t.fLow = t.minimum_frequency_hz;
t.fHigh = t.maximum_frequency_hz;

% LFDCS's min/max_frequency_hz columns may actually encode sweep
% start/end frequency rather than true numeric min/max -- a downsweep
% D-call would then have "minimum" (the start, higher) > "maximum" (the
% end, lower). Swap so fLow is always numerically <= fHigh, which is
% what downstream QC/SNR code assumes regardless of call direction.
swapIx = t.fLow > t.fHigh;
if any(swapIx)
    tmp = t.fLow(swapIx);
    t.fLow(swapIx) = t.fHigh(swapIx);
    t.fHigh(swapIx) = tmp;
end

t.freq = [t.fLow t.fHigh];

%% File information
if any(strcmpi(t.Properties.VariableNames,'filename'))
    t.BeginFile = string(t.filename);
else
    t.BeginFile = strings(nDetect,1);
end

%% Metadata
% NOTE: wrap in cell *before* repmat (not repmat-then-cellstr) -- repmat
% of an empty string stays empty regardless of nDetect, which breaks the
% table-row-count assignment whenever soundFolder/siteCode/classification
% is '' (the default when not supplied by the caller). Wrapping in {...}
% first replicates a 1x1 cell correctly regardless of the contained
% value's emptiness. Matches the safe pattern already used in
% cnnTableToDetection.m.
t.soundFolder = repmat({soundFolder},nDetect,1);
t.siteCode = repmat({siteCode},nDetect,1);
t.classification = repmat({classification},nDetect,1);

%% Sort detections chronologically
t = sortrows(t,'t0');

end
