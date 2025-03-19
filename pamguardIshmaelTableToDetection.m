function [tab, tableName, tableNum] = ...
    pamguardIshmaelToDetection(...
    databaseFile, databaseTable, soundFolder, siteCode, classification, duration)
%a = pamguardIshSpectroCorrTableToDetection(databaseFile,soundFolder,
%       siteCode,classification, duration)
% Load a Pamguard database table containing Ishmael Spectrogram Correlation
% detections. Convert this into a Matlab struct, 'a', with some basic
% metadata for an acoustic detection. The metadata in Raven Selection
% Tables is highly variable, and here we assume that each Selection Table
% contains the following columns/fields:
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

detectorName = databaseTable;


tableName = showPamguardTables(databaseFile);
% Name of table (all detectors have the same name + a number)

% tableIndex = strncmpi(tableName,detectorName,length(detectorName));
% tableIndex = ~cellfun(@isempty,regexp(tableName,[detectorName '_\d*_']));
tableIndex = contains(tableName,detectorName);
tableName = tableName(tableIndex);

%  Ishmael_spectrogram_correlation_X_events
for i = 1:length(tableName)
    tab{i} = loadPamguardDetectionTable(databaseFile,tableName{i});

    % Add on sound folder, site code & classification, and fix duration
    nDetect = height(tab{i});
    if nDetect==0 % Table was empty
        continue
    end
    tab{i}.soundFolder =  cellstr(repmat(soundFolder,nDetect,1));             % Location of audio files for this detection
    tab{i}.siteCode = cellstr(repmat(siteCode,nDetect,1));                    % Append site to data structure
    tab{i}.classification = repmat({classification},nDetect,1);      % Append classification to data structure;

    % For Pamguard's Ishmael spectrogram correlation detector, the duration
    % is actually that of the kernel, not the amount of time above
    % threshold
    tab{i}.duration = repmat(duration,nDetect,1);

    % Also, the start time should be the peak time minus the duration,
    % rather than first slice above threshold
    tab{i}.t0=tab{i}.t0+(tab{i}.PeakDelaySeconds - duration)/86400;
    tab{i}.tEnd = tab{i}.t0+duration/86400;
end

if nargout == 3
    % Extract numbers (e.g. threshold from Pamguard database table names %%
    % Assume the threshold for each detector is embedded in the name of the database
    % tables in Pamguard. For example The table called Ishmael_energy_sum_00001_events
    % should be an Ishmael energy detector module with a threshold of 1.
    for i = 1:length(tableName)
        try
            [intPart, decPart ] = strread(tableName{i},[detectorName '_%d_%d_events']);
            tableNum(i) = str2num([num2str(intPart) '.' num2str(decPart)]);
        catch
            try
                threshString{i} = strread(tableName{i},[detectorName '_%d_events']);
                tableNum(i) = threshString{i};
            catch
                threshString{i} = '';
                tableNum(i) = nan;
            end
        end
    end
end




