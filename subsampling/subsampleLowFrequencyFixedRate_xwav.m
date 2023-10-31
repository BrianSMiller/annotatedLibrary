function subsampleLowFrequencyFixedRate(fileInfo,params);
if nargin<1  
    %%   
    params = createExampleParams
    example = evalc(['disp(params)']);
    error('The first argument for this function must be a valid parameter structure.\nFor example:\n%s',example);
end
    
%% No user-adjustable parameters below this line. Edit below at your own discretion
params.fileInfo = fileInfo;

if params.constantSubSampleRate % distribute chunks over the 353-day duration of the Kerguelen 2014 dataset
    numDaysPerSite = 353;     
else % distribute chunks over the actual duration of the recording
    numDaysPerSite = params.endDate - params.startDate; 
end

numHoursPerSite = numDaysPerSite * 24;
sampleSpacing = numDaysPerSite/params.numberOfChunks * 24; % sample spacing in hours

startTime = rem(params.startHour,sampleSpacing)/24; % To be added to 1st day of sampling

%% Generate vector of evenly spaced samples for this site
subsampleStartTime = nan(params.numberOfChunks,1);

% Calculate start and end dates
startDate = params.startDate+startTime;
endDate = [params.startDate] + min(numDaysPerSite);
temp = datevec(startDate:sampleSpacing/24:endDate);

% Set minutes and seconds to zero
temp(:,[5,6]) = zeros(length(temp),2);
subsampleStartTime(:,1) = datenum(temp);
subsampleEndTime(:,1) = subsampleStartTime + params.durationOfChunk/86400;

%% 
% Check whether there is any temporal aliasing or missing hours,days,months
[year, month, day, hour] = datevec(subsampleStartTime);
figure;
subplot(3,1,1);
hist(month,[1:12]);
xlabel('month');
axis tight;

subplot(3,1,2);
hist(day,[1:31]);
xlabel('day');
ylabel('frequency of occurrence');
axis tight;

subplot(3,1,3);
hist(hour,[0:23]);
xlabel('hour');
axis tight;
%%
% Check the relationship between chosen days, hours, and months.  The most
% important plot here is the first one of hours against day of the month.
% This helps to flag if there is any temporal repetetive patterning to be
% aware of.
figure;
subplot(3,1,1);
plot(day,hour,'o');
xlabel('day of month');
ylabel('hours');
axis tight;

subplot(3,1,3);
plot(month,day,'o');
xlabel('month');
ylabel('day');
axis tight;

subplot(3,1,2);
plot(month,hour,'o');
xlabel('month');
ylabel('hours');
axis tight;
%% Extract sub-samples from data and print the metadata to the console

% Make sure the output folder exists, and create a new subfolder with the
% code for this recording
[s, mess, messid] = mkdir(params.outputFolder, params.code);
if s == 0
    warning('Output folder does not exist. Aborting');
    if params.createOutputWavFiles
        return
    end
end

% The subsample can be divided into smaller portions that each span the 
% entire dataset, but with less temporal resolution. This can be useful if
% might only have enough time to analyse 100 files instead of 200.
% Any Nth subsample can be chosen as A:N:params.numberOfChunks, where A=1,2,...,N
% For example to divide the subsample into two:
% First half of subsamples chosen with the following for statement:
%     for j = 1:2:params.numberOfChunks
% Second half of chosen as 
%     for j = 2:2:params.numberOfChunks
% Otherwise all subsamples will be selected with the following statement:
%     for j = 1:1:params.numberOfChunks
fid = fopen([params.code '_subsample.csv'],'a');
for j = 1:1:params.numberOfChunks;
    timeString = datestr(subsampleStartTime(j),'yyyymmdd_HHMMSS');
    newFilename = [params.outputFolder filesep params.code filesep timeString '.wav'];
    
    % Write a new output file if it doesn't already exist
    if params.createOutputWavFiles;
        % Check that output folder exists, and create it if it doesn't
        if (~exist(newFilename,'file'))
            [sound, ~, wavInfo] = getAudioFromFiles(params.fileInfo,...
                subsampleStartTime(j),subsampleEndTime(j));
%             [sound, sampleRate] = audioread(data(i).fileInfo(fileIndex).fname);
            if isempty(wavInfo)
                warning('Warning: No audio found between %s and %s',...
                    datestr(subsampleStartTime(j)),datestr(subsampleEndTime(j)));
                continue;
            end
            newSound = decimate(sound,wavInfo(1).sampleRate/params.outputSampleRate);
            audiowrite(newFilename,newSound,params.outputSampleRate);
        end
    end
    
    % Print a table of the subsample metadata to the console. Output can be
    % pasted into an Excel spreadsheet to generate a metadata file.
    fprintf('%g\t%s\t%s\t%s\t%s\n',...
        j,...
        params.code,...
        newFilename,...
        datestr(subsampleStartTime(j,1),'yyyy-mm-dd HH:MM:SS'),...
        datestr(subsampleEndTime(j,1),'yyyy-mm-dd HH:MM:SS')...
        );
    
    % Print a table of the subsample metadata to a CSV file to assist in
    % tracking which files have been manually annotated

        fprintf(fid,'%g,%s,%s,%s,%s\n',...
        j,...
        params.code,...
        newFilename,...
        datestr(subsampleStartTime(j,1),'yyyy-mm-dd HH:MM:SS'),...
        datestr(subsampleEndTime(j,1),'yyyy-mm-dd HH:MM:SS')...
        );
end
fclose(fid);

function params = createExampleParams
% data.code will be used when creating output files and folders
    params.code = 'casey2014';
    
    % data.inputFolder should be edited to match the location of your wav files
    params.inputFolder = 'M:\Casey2014\';
    
    params.timeStampFormat = 'yyyy-mm-dd_HH-MM-SS';
    
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