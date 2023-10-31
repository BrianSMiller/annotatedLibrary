%% Generate a list of hours from Casey2014 and Kerguelen2014
% Sample data using a fixed sample rate with a random starting time 

% 10% of the analysis should be blind replicates drawn from already
% analysed data in order to obtain a measure of reproducibility/variability
% of the analyst. These blind replicates will need to be carefully selected 
% from data that had already been analysed so that it will contain a
% statistically relevant number of the calls of interest (e.g >40), but
% also will have some periods with no or few calls.

sitecode = {'casey2014','kerguelen2014'};
refreshFileInfo = false;
copyWavFiles = false;
totalNumberOfChunks = 2000; % hour long chunks 
blindResamples = .10; % Reserve this percentage of hours for blind
                      % resampling to investigate variability of analyst 

newWavSampleRate = 1000; % Resample wav files to this rate (Hz) when copying them
newWavFolder = '\\aad.gov.au\files\ftproot\Public\BrianMiller\sorp\'; % location to copy new files


%% No changs to file required below this line
for i = 1:length(sitecode);
    data(i) = loadRecorderMetaData(sitecode{i});
    if refreshFileInfo
        wavFolderInfo(data(i).wavFolder);
        save([data(i).code 'fileInfo.mat'],'fileInfo');
    end
end

for i = 1:length(sitecode);
    load([data(i).code 'fileInfo.mat'],'fileInfo');
    data(i).fileInfo = fileInfo;
end

numDaysPerSite = [data.endDate] - [data.startDate];
numDaysPerSite = 365;
numHoursPerSite = numDaysPerSite * 24;
numSites = length(data);    % total numbe of sites

totalSamplesPerSite = totalNumberOfChunks/numSites; 
numSamplesPerSite = totalSamplesPerSite - totalSamplesPerSite*blindResamples;
sampleSpacing = min(numDaysPerSite)/numSamplesPerSite * 24; % sample rate in hours

% Pick a random starting point
startHour = randsample(floor(min(numHoursPerSite)),length(data)); 
startHour = [6916, 816];% Output from above, but hard-coded 
                        % here to be able to regenerate the spreadsheet below
startTime = rem(startHour,sampleSpacing)/24; % To be added to 1st day of sampling 

%% Generate vector of evenly spaced samples for each site
sampleStartTime = nan(numSamplesPerSite,numSites);
for i = 1:length(data);
    % Calculate start and end dates
    startDate = data(i).startDate+startTime;
    endDate = [data(i).startDate] + min(numDaysPerSite);
    temp = datevec(startDate:sampleSpacing/24:endDate);
    
    % Eliminate minutes and seconds 
    temp(:,[5,6]) = zeros(length(temp),numSites);
    sampleStartTime(:,i) = datenum(temp);
%     sampleStartTime(:,i) = dateround(sampleStartTime(:,i),'hour','round')
end

%% Within each day sample randomly sample X hours

% Initial 100 files chosen as j = 1:9:numSamplesPerSite
% Second 100 files chosen as j = 5:9:numSamplesPerSite
for i = 1:numSites
    for j = 1:9:numSamplesPerSite;

           fileTimeDiff = abs([data(i).fileInfo.startDate]-sampleStartTime(j,i));
           [timeDiff, fileIndex] = min(fileTimeDiff);
           if timeDiff > 2/24;
               continue;
           end
           [pathname, filename, ext] = fileparts(data(i).fileInfo(fileIndex).fname);
           newFilename = [newWavFolder filesep data(i).code filesep filename '.wav'];
           
           % Write a new downsample file if it doesn't already exist
           if copyWavFiles;
           if (~exist(newFilename,'file'))
               [sound, sampleRate] = audioread(data(i).fileInfo(fileIndex).fname);
               newSound = decimate(sound,sampleRate/newWavSampleRate);
               audiowrite(newFilename,newSound,newWavSampleRate);
           end
        end
        
        fprintf('%s\t%s\t%s\t%s\n',...
            data(i).code,...
            filename,...
            datestr(sampleStartTime(j,i),'yyyy-mm-dd HH:MM:SS'),...
            datestr(sampleStartTime(j,i) + 1/24,'yyyy-mm-dd HH:MM:SS')...
        );
    end    
end

%% Locate wav files, downsample, and copy to a new location for download




