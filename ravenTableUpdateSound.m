function d = ravenTableUpdateSound(dataFile, newSoundInfo, newFileName)
% Update a Raven Selection table that was created from a different sound so
% that annotations display correctly with the new sound. E.g. Raven
% annotations were created using a subset of the annotated data, but we
% want to view them within the full dataset (or vice-versa). Simply
% importing the selection table from the subset into the Raven "sound" that
% contains the full dataset will not work correctly. This is because Raven
% calculates the time of the Annotation relative to the duration of all the
% files within the original Raven "sound", rather than by using the
% absolute time or the number of seconds within each actual audio file. To
% put it another way, Raven treats a folder of sound files as if it were a
% single large continuous sound. 
% 
% NB: Raven table must have BeginDateTime for this script to work 
%
% Here we use the function wavFolderInfo from the SoundFolder package to
% create a new modified selection table that will work in the new Raven
% "sound"

if nargin < 3
    [outpath outname outext] = fileparts(dataFile);
    newFileName = [fullfile(outpath,[outname '_updated']), outext];
end

d = readtable(dataFile,'delimiter','\t');

t = datenum(d.BeginDateTime); 
% Remove observations that are Out of bounds
oob = find(t<min([newSoundInfo.startDate]) | t>max([newSoundInfo.endDate]));
d(oob,:)=[];
t(oob)=[];

% Lookup the file index for each detection using histc
[~, index] = histc(t,[newSoundInfo.startDate]');

% No valid detections in this file
if isempty(index)
    copyfile(dataFile,newFileName);
    return;
end
oob = (index < 1);
index(oob) = [];
d(oob,:) = [];
t(oob) = [];

% First we compute the amount of seconds since the start of the file
secSinceFileStart = (t-[newSoundInfo(index).startDate]')*86400;

% Assume all are same duration
oob = find(secSinceFileStart>newSoundInfo(1).duration);
d(oob,:)=[];
t(oob)=[];
secSinceFileStart(oob) = [];    

% Lookup the file index for each detection using histc
[~, index] = histc(t,[newSoundInfo.startDate]');


% Now add in the cumulative number of seconds since the start of the sound
newFileEndSec = cumsum([newSoundInfo.duration]);
newFileStartSec = [0 newFileEndSec(1:end-1)];
d.BeginTime_s_ = secSinceFileStart + newFileStartSec(index)';
d.EndTime_s_ = (d.BeginTime_s_ + d.DeltaTime_s_); 

% Convert seconds since "sound" start, to samples 
sampleRate = [newSoundInfo(index).sampleRate]';
d.BegFileSamp_samples_ = floor(secSinceFileStart .* sampleRate);
endSecSinceFileStart = secSinceFileStart + d.DeltaTime_s_;
d.EndFileSamp_samples_ = floor(endSecSinceFileStart .* sampleRate);

% Raven does not want the full path in it's file names, so remove this
BeginFile = {newSoundInfo(index).fname}';
[~,fileName,ext] = cellfun(@fileparts,BeginFile,'UniformOutput',false);
d.BeginFile = strcat(fileName,ext);

[~, endIx] = histc(t+(d.DeltaTime_s_/86400),[newSoundInfo.startDate]');
EndFile = {newSoundInfo(endIx).fname}';
[~,fileName,ext] = cellfun(@fileparts,BeginFile,'UniformOutput',false);
d.EndFile = strcat(fileName,ext);

% Save the modified selection table 
fid = fopen(dataFile,'r');
header = fgetl(fid);
fclose(fid);

fid = fopen('RavenTableHeader.txt','w');
fprintf(fid,'%s\n',header);
fclose(fid);

% Matlab has changed the names of the header, so we discard those
writetable(d,'RavenTableBody.txt','delimiter','\t','WriteVariableNames',false);
cmd = sprintf('!type RavenTableHeader.txt RavenTableBody.txt > "%s"',newFileName);
eval(cmd);
