%% Subsample data from Casey2014 long term recording 
% Sample data using a fixed sample rate with a random starting time. The
% idea is to obtain a representative subsample across time and noise
% conditions.

% data.code will be used when creating output files and folders
params.code = 'kerguelen2014'; 

% data.inputFolder should be edited to match the location of your wav files
params.inputFolder = 'M:\kerguelen2014\'; 

% Output files will go into this location
params.outputFolder = '\\aad.gov.au\files\ftproot\Public\BrianMiller\sorp\test\'; 

% Total number of discrete subsampling periods (i.e. chunks). 
params.numberOfChunks = 200;

% The duration of each chunk in seconds
params.durationOfChunk = 3600;  

% Resample wav files to this rate (Hz) if creating output files
params.outputSampleRate = 1000; 

% If refreshFileInfo is true, then the wav metadata and timestamps will
% always be loaded from disk. Otherwise, the wav metadata and timestamps
% will be cached so that they load quickly.
params.refreshFileInfo = false;

% If createOutputWavFiles is true, new wav files will be created, otherwise
% the filenames, start and stop times are printed to the console
params.createOutputWavFiles = true; 

% If constantRate is true, then the subsample rate will be the number of
% numberOfChunks/(365*24). This may result in fewer than numberOfChunks sub
% samples if the duration of recording is less than 1 full year.
% If constantRate is false, the total number of chunks will be evenly 
% distributed throughout the full duration of the recording 
params.constantSubSampleRate = true;

% Timestamp format can be used  by wavFolderInfo
params.timeStampFormat = 'yyyy-mm-dd_HH-MM-SS';

% Pick a random starting point
startHour = floor(rand(1)*25)

% Output from above, but hard-coded here to be able to reproduce the result
params.startHour = 816; 

% Load fileinfo and pass to subsample script allows datasets other than
% wavFiles
fileInfo = wavFolderInfo(params.inputFolder, params.timeStampFormat, params.refreshFileInfo);
params.startDate = min([fileInfo.startDate]); 
params.endDate = max([fileInfo.endDate]);
