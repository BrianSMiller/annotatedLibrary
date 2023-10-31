function [tab, tableName] = ...
    pamguardIshSpectroCorrTableToDetection(...
        databaseFile, soundFolder, siteCode, classification, duration)
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

detectorName = 'Ishmael_spectrogram_correlation'; 


tableName = showPamguardTables(databaseFile);
% Name of table (all detectors have the same name + a number) 

tableIndex = strncmpi(tableName,detectorName,length(detectorName));
tableName = tableName(tableIndex);

%  Ishmael_spectrogram_correlation_X_events
for i = 1:length(tableName);
    tab{i} = loadPamguardDetectionTable(databaseFile,tableName{i});
    
    % Add on sound folder, site code & classification, and fix duration
    nDetect = height(tab{i});
    
    tab{i}.soundFolder =  cellstr(repmat(soundFolder,nDetect,1));             % Location of audio files for this detection
    tab{i}.siteCode = cellstr(repmat(siteCode,nDetect,1));                    % Append site to data structure
    tab{i}.classification = repmat({classification},nDetect,1);      % Append classification to data structure;

    % For Pamguard's Ishmael spectrogram correlation detector, the duration
    % is actually that of the kernel, not the amount of time above
    % threshold
    tab{i}.duration = repmat(duration,nDetect,1);
    
    % Also, the start time should be the peak time minus the duration,
    % rather than first slice above threshold
    tab{i}.t0=tab{i}.t0+tab{i}.PeakDelaySeconds/86400;
    tab{i}.tEnd = tab{i}.t0+duration/86400;    
end




