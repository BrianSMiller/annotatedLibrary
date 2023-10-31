%% Generate a list of hours from Casey2014 and Kerguelen2014
% Sample data using a fixed daily sample rate and randomly sampling 9 hours
% each day
% Use a fixed rate sampling scheme with a random starting time to we would
% investigate 9 randomly sampled hours per day approximately every 3.5 days. 
%
% 10% of the analysis should be blind replicates drawn from already
% analysed data in order to obtain a measure of reproducibility/variability
% of the analyst. These blind replicates will need to be carefully selected 
% from data that had already been analysed so that it will contain a
% statistically relevant number of the calls of interest (e.g >40), but
% also will have some periods with no or few calls.

sitecode = {'casey2014','kerguelen2014'};
totalNumberOfChunks = 2000; % hour long chunks 
numDaysToSample = 100;
blindResamples = .10; % Reserve this percentage of hours for blind
                      % resampling to investigate variability of analyst 

for i = 1:length(sitecode);
    data(i) = loadRecorderMetaData(sitecode{i});
end

numDaysPerSite = floor([data.endDate] - [data.startDate]);
sampleSpacing = min(numDaysPerSite)/numDaysToSample; % 1st tier sample rate in days

numSites = length(data);    % total numbe of sites
totalSamplesPerSite = totalNumberOfChunks/numSites; 
samplesPerSite = totalSamplesPerSite - totalSamplesPerSite*blindResamples;
samplesPerDay = samplesPerSite/numDaysToSample;

% Pick a random starting point
startDay = randsample(350,length(data)); 
startDay = [336,56]; % initially chosen randomly from above, but hard-coded 
                     % here to be able to regenerate the spreadsheet below
startDay = rem(startDay,sampleSpacing); % To be added to 1st day of sampling 

%% Generate a vector of evenly spaced days to sample
daysSampled = nan(numDaysToSample,numSites);
for i = 1:length(data);
    endDate = [data(i).startDate] + min(numDaysPerSite);
    temp = datevec(data(i).startDate+startDay:sampleSpacing:endDate);
    temp(:,[5,6]) = zeros(length(temp),2); % zero out minutes and seconds
    daysSampled(:,i) = datenum(temp);
end

%% Within each day sample randomly sample X hours
for j = 1:length(daysSampled);
%     while
    fprintf('%s, %s, %s,\n',...
        data.code{i},...
        datestr(data.startTime(j),'yyyy-mm-dd HH:MM:SS'),...
        datestr(stopTime(j),'yyyy-mm-dd HH:MM:SS')...
        );
%     end
end

