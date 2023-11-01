function [ tab ] = loadPamguardDetectionTable( fileName, tableName, verbose );
% loadPamguardDetectionTimes: Load detection times of clicks from PAMGuard's 
% database.
% fileName - name of PAMGuard's sqlite database
% tableName - name of the database table that contains click detections

% folder = 'Z:\Acoustics\workInProgress\kerguelen2014_spermWhales\';
% database = 'kerguelen2014_spermWhale.sqlite3';
% fileName = [folder database];
if nargin < 2
    tableName = 'Click_Detector_Clicks';
    fprintf('No database table specified. Using default table: %s\n',tableName);
end

if nargin < 3
    verbose = false;
end

julianDayOffset = 1.721058500000000e+06; % Number of julian days before (0000 Jan 00)

dbid = mksqlite('open',fileName);
tables = mksqlite('show tables');
if verbose
    fprintf('Opened sqlite database: %s\nLoading detection times from table: %s.\n',fileName, tableName);
    fprintf('Still working. Loading large databases can take a few minutes...\n');
    tic;
end
result = sql(dbid,'select * from %s order by UTC ',tableName);
mksqlite(dbid,'close');
if verbose
    fprintf('Finished reading from database.\n');
    toc
end
tab = struct2table(result);
try
tab.t0 = datenum(tab.UTC,'yyyy-mm-dd HH:MM:SS.fff');
tab.tEnd = tab.t0 + tab.DurationSeconds/86400;
tab.fLow = tab.lowFreq;
tab.fHigh= tab.highFreq;
tab.freq = [[tab.fLow] [tab.fHigh]];  % Frequency vector in Hz
catch
end
    

% fprintf('Converting timestamps to Matlab datenum.\n');
% fprintf('This could also take a while...\n');
% tic; detectionTimes = datenum({result.UTC},'yyyy-mm-dd HH:MM:SS.fff'); toc
% detectionTimes = [result.julianday_UTC_] - julianDayOffset;


