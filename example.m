%% Example analysis for annotatedLibrary
% The following example illustrates basic usage of some of the commonly 
% used functions from the annotatedLibrary github repository 
% https://github.com/BrianSMiller/annotatedLibrary.


%% Define local locations of data (wav files and annotations)

sorpFolder = 'w:\annotatedLibrary\SORP\';
siteCode = 'kerguelen2014';
ravenFile = [sorpFolder, siteCode, '\Kerguelen2014.Bm.Ant-Z.selections.txt'];

%% Register a new wavFolder and load metadata

wavFolder = [sorpFolder siteCode '\wav\'];
[~, timestampFormat] = guessFileNameTimestamp(wavFolder);
wavInfo = wavFolderInfo(wavFolder,timestampFormat);

%% Load a Raven Selection Table

% Simple wrapper around readtable to link the selection & soundFolder
rt = ravenTableToDetection(ravenFile,wavFolder,siteCode,'Z-call');


%% Estimate SNR for annotations 
params.showClips = false;
params.noiseDelay = 1;
params.freq = [24 29];  % Antarctic blue whale: Unit-A of full Z-calls

[snr, rt.signalRms, rt.noiseRms, rt.noiseVar] = annotationSNR(rt,params);
rt.snr = 10*log10(rt.signalRms) - 10*log10(rt.noiseRms);

% Plot results as a time series
figure;
scatter(rt.t0,rt.snr);
datetick('x');
ylabel('SNR (dB)');
grid on;


%% Simple 'sensitivity analysis' to see if a parameter is sensitive to SNR

thresholds = -6:3:14;

tiledlayout(2,1);
ax2 = nexttile;
histogram(rt.snr,thresholds);
xlabel('SNR (dB)')
ylabel('Number of detections')
grid on;

% Calculate mean, std, 25, 50, 75 percentiles of durations for each
% threshold
duration = [];

for i = 1:length(thresholds);
    ix = rt.snr > thresholds(i);
    duration.Q1(i) = quantile(rt.Dur90__s_(ix),0.25);
    duration.median(i) = median(rt.Dur90__s_(ix));
    duration.Q3(i) = quantile(rt.Dur90__s_(ix),0.75);
    duration.mean(i) = mean(rt.Dur90__s_(ix));
    duration.std(i) = std(rt.Dur90__s_(ix));
end

ax3 = nexttile;
errorbar(thresholds,duration.mean,duration.std);
hold on;
plot(thresholds,duration.median,'ok');
line([thresholds;thresholds],[duration.Q1; duration.Q3], ...
    'lineWidth',2,'color','k');
xlabel('SNR (dB)');
ylabel('Duration (s) [90% energy]')
grid on;
linkaxes([ax2 ax3],'x');

%% Sensitivity analysis conclusions 
% On the surface, the results of the sensitivity analysis are a bit
% ambiguous. It looks to me like the duration appears to be increasing with
% SNR at 3 dB as well as at 12 dB. 
% 
% However, there are few annotations with SNR > 12 dB, so I wouldn't rely
% on that last point for anything. And the increase at threshold of 3 dB is
% only from 10.37 - 10.78 = 0.41 s on top of a duration of about 10 s (so
% around 4%). Since the effect appears small, but potentially could be
% real,statistical analysis might be warranted to confirm or reject any
% relationship between SNR and duration. 
% 
% Additionally, it's worth noting that the outcome here could also be due
% to the nature of the dataset. For example, there may be a relationship
% between SNR and duration that we are unable to detect because we have
% only a limited range of SNRs from -6 to 13.3 dB. The relationship may be
% more apparent if we had substantially more low- and/or high-SNR 
% detections.
%


