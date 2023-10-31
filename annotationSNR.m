function [snr, rmsSignal, rmsNoise, noiseVar, fileInfo] = annotationSNR(annot,params)
% Measure the signal-to-noise ratio (SNR) of a detection from a Raven
% selection table.
%
% SNR is defined here following equation 6.26 in Lurton (2010) "An
% Introduction to Underwater Acoustics, Principles and Applications. 2nd
% Edition". The noise floor is computed from a segment of time just before
% and after the event. The time segment containing noise will have the same
% bandwidth and duration as the event.
%
% ANNOTATION is a data structure that contains a start and end time and
% duration (t0, tEnd), and duration in and frequency vector (freq) of high
% and low frequencies. t0 and tEnd must be Matlab datenums. Frequency
% vector should be in Hz. The detection must also have a field called
% soundFolder, that is the path to a folder of wav files (i.e. can be read
% via the function wavFolderInfo.m).
%
% PARAMS is a data structure with additional fields that can influence the
% way that the SNR is measured. At the moment.
%   PARAMS.noiseDelay - is a scalar representing the time gap, in seconds,
%   between the noise and signal measurements.
%   PARAMS.freq - is a 2x1 vector defining the lower and upper frequency
%   band to use for SNR measurements. If params.freq is not specified, then
%   the frequency band from ANNOTATION.freq will be used for SNR
%   measurements
%   PARAMS.showClips will display spectrogram with marks around noise and
%   signal.
%   PARAMS.spectroParams is a sub-structure that can alter the appearance
%   of the spectrogram (only if PARAMS.showClips is true). By default,
%   spectrogram parameters are geared towards low frequency sounds 1 s
%   long time slices with a lot of overlap.
%     spectroParams.pre = 4; % pre-detection sound buffer (seconds)
%     spectroParams.post = 4;% post-detection sound buffer (seconds)
%     spectroParams.win = sampleRate; % FFT & window length (samples)
%     spectroParams.overlap = floor(sampleRate*0.85); % Time slice overlap
%     (samples)
%     spectroParams.yLims = [0 125]; % Spectrogram frequency limits (Hz)
%
%% Brian Miller, Australian Antarctic Division 2017.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process params and set some default parameters if any required params are
% missing.
% Toggle showClips to plot spectrograms of each annotation * noise
showClips = false;

if nargin < 2
    params.noiseDelay = 0;
    params.showClips = false;
    params.freq = [0 500];
end

if isfield(params,'freq')
    freq = params.freq;
else
    % User has specified freq in the function parameters
    if isstruct(annot) && length(annot)==1
        % Use the frequency limits of the annotation
        freq = annot.freq;
    end
end

if isfield(params,'showClips')
    showClips = params.showClips;
end

if ~isfield(params,'noiseDelay')
    params.noiseDelay = 0;
end

if ~isfield(params,'noiseDuration')
    % noiseDuration = 'randomBeforeAndAfter';
    % noiseDuration = '1minuteBeforeAndAfter';
    % noiseDuration = '30sBeforeAndAfter';
    params.noiseDuration = 'beforeAndAfter';
end

if ~isfield(params,'snrType')
    params.snrType = 'spectrogram';
    % params.snrType = 'quantiles';
    % params.snrType = 'timeDomain';
end

if ~isfield(params,'metadata')
    params.metadata = [];
end

% TODO make these spectrogram parameters user-adjustable input-parameters
% desired number of spectrogram slices (actual number will be slightly
% different depending on overlap and duration of clip)
nSlices = 30; 
% proportion between (0, 1] of overlap between spectrogram slices. 0 means
% no overlap between slices, and 1 would mean full overlap.
overlap = 0.75; 


%%%% End of processing of input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Main body of function below

% If annotation is a table, then convert it to a struct before proceeding
if istable(annot)
    annot = table2struct(annot);
end

% If annotation is a vector (of many annotations) loop through them all
if length(annot) > 1

    snr = nan(size(annot));
    rmsSignal = snr; rmsNoise = snr; noiseVar = snr; str = ''; tic;
    for i = 1:length(annot)
        [snr(i), rmsSignal(i), rmsNoise(i), noiseVar(i)] = annotationSNR(annot(i),params);
        if rem(i,100)==0 || i==length(annot)
            fprintf(repmat('\b',1,length(str)));
            str = sprintf('%g/%g SNR measurements completed in %4.1f s\n',...
                i,length(annot),toc);

            fprintf('%s',str);
        end
    end
    return
end

soundFolder = wavFolderInfo(annot.soundFolder);
[annot.audio, ~, annot.fileInfo] = getAudioFromFiles(soundFolder,annot.t0,annot.tEnd);

fileInfo = annot.fileInfo;
if isempty(annot.fileInfo) || isempty(annot.audio)
    snr = nan;
    rmsSignal = nan;
    rmsNoise = nan;
    noiseVar = nan;
    return;
end

% Assume all files are at the same sample rate
sampleRate = annot.fileInfo(1).sampleRate;
nfft = 2^(nextpow2(floor(annot.duration/nSlices / overlap * sampleRate)));
nOverlap = floor(nfft*overlap); % overlap in samples


% Create a noise data structure just a bit earlier than the annotation
noise = annot;
switch(params.noiseDuration)
    case 'before'
        noise.tEnd = annot.t0-params.noiseDelay/86400; % end of noise is 1 s before annotation
        noise.t0 = noise.tEnd-annot.duration/86400;
        [noise.audio, ~, noise.fileInfo] = getAudioFromFiles(soundFolder,noise.t0,noise.tEnd);
        excludeTimes = [];
    case '25sBefore' % as requested by Franciele Castro for SORP ATWG post-doc
        % load only 25 s before the detection
        noiseDuration = 25/86400; % 25 s (converted to days/matlab datenum)
        noise.tEnd = annot.t0-params.noiseDelay/86400; % end of noise is 1 s before annotation
        noise.t0 = noise.tEnd - noiseDuration;
        excludeTimes = []; % Placeholder for logic to check for other signals
        [noise.audio, ~, noise.fileInfo] = getAudioFromFiles(soundFolder,noise.t0,noise.tEnd);
    case '30sBeforeAndAfter'
        noise.tEnd = annot.tEnd + (30 + params.noiseDelay)/86400;
        noise.t0 = annot.t0 - (30 + params.noiseDelay)/86400;
        excludeTimes = [annot.t0 annot.tEnd] + params.noiseDelay*([-1 1]/86400);
        [noise.audio, ~, noise.fileInfo] = getAudioFromFiles(soundFolder,noise.t0,noise.tEnd,excludeTimes);
    case 'randomBeforeAndAfter'
        params.noiseDelay = rand * params.noiseDelay;
        noise.tEnd = annot.tEnd + (0.5*annot.duration + params.noiseDelay)/86400;
        noise.t0 = annot.t0 - (0.5*annot.duration + params.noiseDelay)/86400;
        excludeTimes = [annot.t0 annot.tEnd] + params.noiseDelay*([-1 1]/86400);
        [noise.audio, ~, noise.fileInfo] = getAudioFromFiles(soundFolder,noise.t0,noise.tEnd,excludeTimes);
    otherwise % Default to using noise before and after the sound
        noise.tEnd = annot.tEnd + (0.5*annot.duration + params.noiseDelay)/86400;
        noise.t0 = annot.t0 - (0.5*annot.duration + params.noiseDelay)/86400;
        excludeTimes = [annot.t0 annot.tEnd] + params.noiseDelay*([-1 1]/86400);
        [noise.audio, ~, noise.fileInfo] = getAudioFromFiles(soundFolder,noise.t0,noise.tEnd,excludeTimes);
end

switch params.snrType
    case 'spectrogram'
        %% Measure signal power and variance using spectrogram *cells*
        [rmsSignal, ~, ~, ~] = spectrogramPowerAndVariance(...
            annot.audio, nfft, nOverlap, nfft, sampleRate, freq,...
            params.metadata);
        [rmsNoise, noiseVar, ~, ~] = spectrogramPowerAndVariance(...
            noise.audio, nfft, nOverlap, nfft, sampleRate, freq,...
            params.metadata);
        annot.rmsLevel = 10*log10(rmsSignal);
        noise.rmsLevel = 10*log10(rmsNoise);

        snr = (rmsSignal - rmsNoise).^2/noiseVar;
        snr = 10*log10(abs(snr));
    case 'spectrogramSlices'
        %% Measure signal power and variance from spectrogram *slices* 
        % This was the method used for the Annotated Library publication
        % (Miller et al 2021), and Deep Learning D-call paper (Miller et al
        % 2022)

        [rmsSignal, ~, ~, ~] = slicePowerAndVariance(...
            annot.audio, nfft, nOverlap, nfft, sampleRate, freq,...
            params.metadata);
        [rmsNoise, noiseVar, ~, ~] = slicePowerAndVariance(...
            noise.audio, nfft, nOverlap, nfft, sampleRate, freq,...
            params.metadata);
        annot.rmsLevel = 10*log10(rmsSignal);
        noise.rmsLevel = 10*log10(rmsNoise);

        snr = (rmsSignal - rmsNoise).^2/noiseVar;
        snr = 10*log10(abs(snr));

    case 'quantiles'
        % Subtract Ln_noise from Ln_signal (e.g. 85th-15th percentile)
        [rmsSignal, rmsNoise, noiseVar] = snrFromQuantile(...
            annot.audio, nfft, nOverlap, nfft, sampleRate, freq);

        annot.rmsLevel = 10*log10(rmsSignal);
        noise.rmsLevel = 10*log10(rmsNoise);

        snr = annot.rmsLevel - noise.rmsLevel;
    case 'timeDomain' %% measure signal power in time domain
        % Here we bandpass filter the time domain signal and noise. We
        % calculate the RMS power of each, and we calculate the variance of
        % the noise and signal to follow eq. 6.25 from Lurton (2010)
        % Introduction to Underwater Acoustics to estimate SNR. 

        % signalPower = bandpower(annotation.audio, sampleRate, freq);
        % noisePower = bandpower(noise.audio, sampleRate, freq);
        cutoffFreq1 = freq(1);
        cutoffFreq2 = freq(2);
        persistent lastCutoffFreq1
        persistent lastCutoffFreq2
        persistent lastSampleRate
        persistent lastd
        %cache filter for improved performance
        if isempty(lastCutoffFreq1) ...
                || isempty(lastCutoffFreq2)...
                || isempty(lastSampleRate) ...
                || cutoffFreq1~=lastCutoffFreq1...
                || cutoffFreq2~=lastCutoffFreq2...
                || lastSampleRate~=sampleRate

            d = designfilt('bandpassfir','FilterOrder',48, ...
                'CutoffFrequency1',40,'CutoffFrequency2',80, ...
                'SampleRate',sampleRate);
        else
            d=lastd;
        end
        lastCutoffFreq1 = cutoffFreq1;
        lastCutoffFreq2 = cutoffFreq2;
        lastSampleRate = sampleRate;
        lastd = d;

        try
            %removeClicks(sourceData, threshold,  power)
            sigFilt = filtfilt(d,annot.audio);
            rmsSignal = rms(sigFilt.^2);

            noiseFilt = filtfilt(d,noise.audio);
            rmsNoise = rms(noiseFilt.^2);
            noiseVar = var(noiseFilt.^2);
        catch me
            sigFilt = nan;
            rmsSignal = nan;
            noiseFilt = nan;
            rmsNoise = nan;
            noiseVar = nan;
            me.stack
        end

        annot.rmsLevel = 10*log10(rmsSignal);
        noise.rmsLevel = 10*log10(rmsNoise);

        % Calculate SNR in dB for a 'power' receiver as per Lurton 2010 eq. 6.25.
        %
        % However this equation for SNR doesn't really make sense if the rmsNoise
        % is greater than rmsSignal, which can happen here since we're using very
        % local estimates of noise.
        %
        % Here we modify Lurton's equation to use the abs(rmsSignal - rmsNoise).
        % This seems a bit wonky, since we're treating the signal as something that
        % supresses noise instead of something that surpasses it, but I can't
        % really think of a better solution...
        snr = (rmsSignal - rmsNoise).^2/noiseVar;
        snr = 10*log10(abs(snr));

end


%%%% Plot the spectrogram of signal & noise. Highlight the signal, noise,
%%%% and SNR.
if showClips
    % Spectrogram parameters included with annotation parameters
    if isfield(params,'spectroParams')
        spectroParams = params.spectroParams;
        if ~isfield(spectroParams,'freq')
            spectroParams.freq = freq;
        end
    else % Create some default spectrogram parameters


        spectroParams.pre = 1;
        spectroParams.post = 1;

        spectroParams.win = floor(sampleRate/4);
        spectroParams.overlapPercent = 0.75;
        spectroParams.overlap = floor(spectroParams.win*spectroParams.overlapPercent);
        spectroParams.yLims = [10 125];
        spectroParams.freq = annot.freq;

        %         binHeight = sampleRate/spectroParams.win;
        binHeight = 1;
        vertScale = 300/diff(spectroParams.yLims);
        spectroParams.verticalPixels =  vertScale * binHeight * diff(spectroParams.yLims);
%         sliceDuration = spectroParams.win/sampleRate*(1-spectroParams.overlapPercent);

        spectroParams.noiseDelay = params.noiseDelay/86400;
        %         pixelsPerSecond = 1;
        %         spectroParams.horizontalPixels = 2*spectroParams.horizontalPixels*c.duration
    end
    horizScale = 1920/860;
    sliceDuration = spectroParams.win/sampleRate*(1-spectroParams.overlapPercent);
    spectroParams.horizontalPixels = horizScale * (3*annot.duration+1+spectroParams.pre+spectroParams.post)/sliceDuration;

    spectroAnnotationAndNoise(annot, noise, soundFolder, spectroParams, snr,...
        params.metadata)
    %%
    if strcmpi(params.snrType,'timeDomain')
        figure(222);
        tN = (0:1/sampleRate:((length(noise.audio)-1)/sampleRate))-length(noise.audio)/sampleRate/2;
        tS = (0:1/sampleRate:((length(annot.audio)-1)/sampleRate));
        plot(tN,noiseFilt.^2,'color',[0.5 0 0]);
        hold on;
        plot(tS,sigFilt.^2,'color',[0 0.5 0]);
        xLim = xlim;
        plot([xLim(1) 0],rmsNoise*[1 1],'color','r','lineWidth',2);
        plot([0 xLim(2)],rmsSignal*[1 1],'color','g','lineWidth',2);
        errorbar(xLim(1),rmsNoise,sqrt(noiseVar),'color',[0 0 0],'lineWidth',3,'linestyle','--');
        %     plot(noiseTime, noisePower,'color',[0.5 0 0]);
        %     hold on;
        %     plot(sigTime+noiseTime(end),sigPower,'color',[0 0.5 0]);
        hold off;
        legend('Noise power','Signal power',...
            'RMS noise','RMS signal','Standard Deviation of Noise power');
        xlabel('Time (s)');
        ylabel('Power (au)');
    end
    pause;
end


%%%% End of Main body of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Internal supporting functions below here

%% Measure power and variance from spectrogram cells
function [power, variance, slicePower, sT] =spectrogramPowerAndVariance(...
    x, window, nOverlap, nfft, sampleRate, freqRange, metadata)
% Measure signal power and variance using short-time fourier tranform
% (spectrogram)
if length(x) < window
    window = length(x);
    nOverlap = 0;
    nfft = window;
end
x = x-mean(x);
[~, sF, sT, specPsd] = spectrogram(x,window,nOverlap,nfft,sampleRate);

% Apply calibration if recording metadata has been supplied
if ~isempty(metadata) 
    caldB = getCalibration(metadata,sF);
    calibration = 10.^(caldB/10);

    specPsd = specPsd./repmat(calibration,size(sT));
end

specPsd = specPsd(sF >= freqRange(1) & sF <= freqRange(2),:);

power = mean(specPsd);
variance = var(specPsd);


% Measure signal power and variance using from slices, rather than cells of
% the spectrogram. This is an attempt to follwow the premise of Lurton 2010 
% of estimating signal power and variance from a time series of 
% band-limited measurements. In particular, we take the mean and variance 
% from the PSDs used to create the spectrogram. For each slice we calculate
% in-band power, and the mean and variance are taken from these. The number
% of samples is the number of spectrogram slices.
function [power, variance, slicePower, sT] = slicePowerAndVariance(...
    x, window, nOverlap, nfft, sampleRate, freqRange, metadata)
% Measure signal power and variance using short-time fourier tranform
% (spectrogram)
if length(x) < window
    window = length(x);
    nOverlap = 0;
    nfft = window;
end
x = x-mean(x);
[~, sF, sT, specPsd] = spectrogram(x,window,nOverlap,nfft,sampleRate);

% Apply calibration if recording metadata has been supplied
if ~isempty(metadata) 
    caldB = getCalibration(metadata,sF);
    calibration = 10.^(caldB/10);

    specPsd = specPsd./repmat(calibration,size(sT));
end

slicePower = nan(size(sT));
for i = 1:length(sT)
    slicePower(i) = bandpower(specPsd(:,i),sF,freqRange,'psd');
end
power = mean(slicePower);
variance = var(slicePower);

function [caldB, freq] = getCalibration(metadata, freq)
% First convert from normalised wav data back to volts
adVpeakdB = 10*log10(1/metadata.adPeakVolt.^2);

% Since we're in frequency domain, correct for the frequency response
% of preamp
frontEndGain_dB = interp1(log10(metadata.frontEndFreq_Hz),...
    metadata.frontEndGain_dB, log10(freq),'linear','extrap');

% Calibration is hydrophone sensitivity + pre-amp + voltage (in dB);
% If using linear units, then we would multiply instead of add
caldB = metadata.hydroSensitivity_dB + frontEndGain_dB - adVpeakdB;

% Fudge the values that are -inf or nan
nanIx = find(isnan(caldB) | isinf(caldB));
caldB(nanIx) = -1000*ones(size(nanIx));


% Measure signal power as 85th percentile of spectrogram intensity
% noise power is 50th percentile of spectrogram intensity
% This method was an experiment to explore a different way of estimating
% SNR and is not recommended for use. 
function [signal, noise, noiseVar] = snrFromQuantile(...
    x, window, nOverlap, nfft, sampleRate, freq)

if length(x) < window
    window = length(x);
    nOverlap = 0;
    nfft = window;
end

[~, F, ~, P] = spectrogram(x,window,nOverlap,nfft,sampleRate);

fIx = F >= freq(1) & F <= freq(2) ;
specPsd = P(fIx,:);
sigIx = specPsd>=quantile(specPsd,0.85,'all');
noiseIx = specPsd<=quantile(specPsd,0.85,'all');
signal = mean(specPsd(sigIx));
noise = mean(specPsd(noiseIx));
noiseVar = var(specPsd(noiseIx));
%         if showClips
%             figure(4);
%             h = pcolor(T,F,10*log10(P)); set(h,'LineStyle','none');
%             hold on;
%             [row, col] = ind2sub(size(specPsd),sigIx);
%             scatter(T(col),fBand(row),'g.');
%             [row, col] = ind2sub(size(specPsd),noiseIx);
%             scatter(T(col),fBand(row),'r.');
%             ylim(freq);
%         figure(4);
%         histogram(specPsd(:));
%         vertline(noise,'r');
%         vertline(signal,'g');
%         end

% Plot the spectrogram showing the detection and the noise time periods
function spectroAnnotationAndNoise(...
    detection, noise, soundFolder, spectroParams, snr, metadata)
for i = 1:length(detection)
    sampleRate = detection.fileInfo(1).sampleRate;
    freq = spectroParams.freq;
    c = detection;
    c.t0 = noise.t0-spectroParams.pre/86400;
    c.tEnd = max([detection.tEnd noise.tEnd]) + spectroParams.post/86400;
    c.duration = (c.tEnd - c.t0)*86400;
    c.audio = getAudioFromFiles(soundFolder,c.t0,c.tEnd);
    if ( c.duration*sampleRate < spectroParams.win ) 
        % Cannot show spectrogram, so just skip
        continue;
    end
    [~, f, t, p] = spectrogram(c.audio,...
        spectroParams.win, spectroParams.overlap, spectroParams.win,...
        sampleRate,'yaxis');
    if ~isempty(metadata)
        caldB = getCalibration(metadata,f);
        calibration = 10.^(caldB/20);
%         plot(20*log10(p./repmat(calibration,size(t)))
        p = p./repmat(calibration,size(t));
    end
    fIx = f >= freq(1) & f <= freq(2);
    pBand = p(fIx,:);
    cLim = quantile(20*log10(pBand(:)),[0 1]);
    figure(111);
    set(gca,'units','pixels');
    %         h = imagesc(t,f, 20 * log10(abs(s)));

    imagesc(t,f, 20*log10(p));
    set(gca,'ydir','normal');
    colormap(flipud(gray))
    set(gca,'clim',cLim);
    set(gca,'xtick',0:floor(c.duration))
    set(gca,'xTickLabel',0:floor(c.duration))
    axPos = get(gca,'position');
    axInset = get(gca,'tightInset');
    pos = get(gcf,'position');
    width = spectroParams.horizontalPixels;
    cb = colorbar('north');
    set(gcf,'pos',[pos(1:2),...
        width+sum(axInset([1,3]))*2.5,...
        spectroParams.verticalPixels+sum(axInset([2,4]))*2]);

    set(gca,'position',[axPos(1:2) width spectroParams.verticalPixels])
    findfigs;

    % Red line for noise
    line(([noise.t0 noise.tEnd]-c.t0)*86400,[1 1]'* freq,...
        'color',[0.5 0 0],'linewidth',2)
    text((noise.t0- c.t0)*86400,freq(1),sprintf('Noise Level = %4.1f dBFS',...
        noise.rmsLevel),'color',[0.5 0 0],'verticalAlignment','top',...
        'BackgroundColor','w');

    % Green line for detections
    line(([detection.t0 detection.tEnd]-c.t0)*86400,[1 1]'* freq,...
        'color',[0 0.5 0],'linewidth',2)
    text((detection.t0 - c.t0)*86400,freq(2),sprintf('Signal Level = %4.1f dBFS',...
        detection.rmsLevel),'color',[0 0.5 0],'verticalAlignment','bottom',...
        'BackgroundColor','w');

    % Black lines for null values (between detection and noise)
    line(([detection.t0 detection.t0-spectroParams.noiseDelay]-c.t0)*86400,...
        [1 1]'*freq,'color','k','linewidth',2);
    line(([detection.tEnd detection.tEnd+spectroParams.noiseDelay]-c.t0)*86400,...
        [1 1]'*freq,'color','k','linewidth',2);

    % SNR
    text((detection.t0- c.t0)*86400,freq(2),sprintf('SNR = %4.1f dBFS    ',snr),...
        'color','g','verticalAlignment','bottom','horizontalAlignment','right');


    ylim(spectroParams.yLims)
    cl = '';
    if isfield(c,'classification'); cl = char(c.classification); end
    title(sprintf('%s - %s',cl, datestr(c.t0) ),'interpreter','none');
    set(gca,'layer','top')
    xlabel('time (s)');
    ylabel('frequency (Hz)');
end
