function [ tables ] = loadPamguardDetectionTable( fileName)
% loadPamguardDetectionTimes: Load detection times of clicks from PAMGuard's 
% database.
% fileName - name of PAMGuard's sqlite database
% tableName - name of the database table that contains click detections

% folder = 'Z:\Acoustics\workInProgress\kerguelen2014_spermWhales\';
% database = 'kerguelen2014_spermWhale.sqlite3';
% fileName = [folder database];

dbid = mksqlite('open',fileName);
result = mksqlite('show tables');
mksqlite(dbid,'close');
tables = {result.tablename};