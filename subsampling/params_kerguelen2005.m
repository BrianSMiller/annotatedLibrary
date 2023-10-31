%% Subsample data from Casey2014 long term recording 
% Sample data using a fixed sample rate with a random starting time. The
% idea is to obtain a representative subsample across time and noise
% conditions.

% data.code will be used when creating output files and folders
params.code = 'kerguelen2005'; 

% data.inputFolder should be edited to match the location of your wav files
params.inputFolder = 'M:\ARPS\Kerguelen05\xwav\'; 

% Output files will go into this location
% params.outputFolder = '\\aad.gov.au\files\ftproot\Public\BrianMiller\sorp\test\'; 
params.outputFolder = 'M:\annotatedLibrary';

% Total number of discrete subsampling periods (i.e. chunks). 
params.numberOfChunks = 200;

% The duration of each chunk in seconds
params.durationOfChunk = 3600;  

% Resample wav files to this rate (Hz) if creating output files
params.outputSampleRate = 500; 

% If refreshFileInfo is true, then the wav metadata and timestamps will
% always be loaded from disk. Otherwise, the wav metadata and timestamps
% will be cached so that they load quickly.
params.refreshFileInfo = false;

% If createOutputWavFiles is true, new wav files will be created, otherwise
% the filenames, start and stop times are printed to the console
params.createOutputWavFiles = false; 

% If constantRate is true, then the subsample rate will be the number of
% numberOfChunks/(365*24). This may result in fewer than numberOfChunks sub
% samples if the duration of recording is less than 1 full year.
% If constantRate is false, the total number of chunks will be evenly 
% distributed throughout the full duration of the recording 
params.constantSubSampleRate = false;

% Timestamp format can be used  by wavFolderInfo
params.timeStampFormat = 'yyyy-mm-dd_HH-MM-SS';

% Pick a random starting point
floor(rand(1,1)*25)
% Output from above, but hard-coded here to be able to reproduce the result
params.startHour = 2;

fileInfo = xwavFolderInfo(params.inputFolder, params.timeStampFormat, params.refreshFileInfo);
% Start and end dates and times of the files as matlab datenums
params.startDate = min([fileInfo.startDate]); % params.startDate = datenum('31-Jan-2005 08:17:06');
params.endDate = max([fileInfo.endDate]); % params.endDate = datenum('25-Feb-2006 02:07:06');
params.startDate = datenum('31-Jan-2005 09:00:00');
params.endDate = params.startDate + 367;