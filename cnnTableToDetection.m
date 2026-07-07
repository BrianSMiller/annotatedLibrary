function t = cnnTableToDetection(ravenFile,soundFolder,...
    siteCode,classification)

% Convert CNN detector table into detection structure
%
% Expected columns:
%
% start_datetime
% end_datetime
% detection_duration_s
% minimum_frequency_hz
% maximum_frequency_hz
% file_start_time
% Probability
% call_type
% filename
%
% Datetime format expected:
% yyyy-MM-dd HH:mm:ss.SSS
%
% Output compatible with downstream callDensity workflow

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

opts = detectImportOptions( ...
    ravenFile,...
    'Delimiter','\t');

% Read datetime columns as text
opts = setvartype( ...
    opts,...
    {'start_datetime','end_datetime','file_start_time'},...
    'char');

t = readtable(ravenFile,opts);

nDetect = height(t);

if nDetect == 0
    return;
end

%% Channel

if any(strcmp(t.Properties.VariableNames,'Channel'))

    t.channel = t.Channel;

else

    t.channel = ones(nDetect,1);

end

%% Datetimes

dt0 = datetime( ...
    t.start_datetime,...
    'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');

dtEnd = datetime( ...
    t.end_datetime,...
    'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');

t.t0   = datenum(dt0);
t.tEnd = datenum(dtEnd);

%% Duration

if any(strcmp(t.Properties.VariableNames,'detection_duration_s'))

    if isnumeric(t.detection_duration_s)

        t.duration = t.detection_duration_s;

    else

        t.duration = str2double(string(t.detection_duration_s));

    end

else

    % Fallback: calculate from timestamps
    t.duration = (t.tEnd - t.t0) * 86400;

end

%% Frequency bounds

if isnumeric(t.minimum_frequency_hz)

    t.fLow  = t.minimum_frequency_hz;
    t.fHigh = t.maximum_frequency_hz;

else

    t.fLow  = str2double(string(t.minimum_frequency_hz));
    t.fHigh = str2double(string(t.maximum_frequency_hz));

end

t.freq = [t.fLow t.fHigh];

%% Optional score field

if any(strcmp(t.Properties.VariableNames,'Probability'))

    if isnumeric(t.Probability)

        t.score = t.Probability;

    else

        t.score = str2double(string(t.Probability));

    end

end

%% Optional call type

if any(strcmp(t.Properties.VariableNames,'call_type'))

    t.callType = string(t.call_type);

end

%% Optional filename

if any(strcmp(t.Properties.VariableNames,'filename'))

    t.BeginFile = string(t.filename);

end

%% Metadata

t.soundFolder = ...
    repmat({soundFolder},nDetect,1);

t.siteCode = ...
    repmat({siteCode},nDetect,1);

t.classification = ...
    repmat({classification},nDetect,1);

%% Sort chronologically

t = sortrows(t,'t0');

end