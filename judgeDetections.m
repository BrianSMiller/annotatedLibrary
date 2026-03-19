function c = judgeDetections(c,sampleRate,nfft,noverlap,pre,post, ...
    startIndex,increment,freqBand,suppressClicks,options)

% JUDGEDETECTIONS  Manually review and label candidate detections
%
%   c = judgeDetections(c, sampleRate, nfft, noverlap) launches an
%   interactive spectrogram viewer for manually labeling detections in a
%   capture history table as true positives or false positives.
%
%   c = judgeDetections(c, sampleRate, nfft, noverlap, pre, post, ...
%       startIndex, increment, freqBand, suppressClicks) uses additional
%       positional options described below.
%
%   c = judgeDetections(___, Name, Value) additionally accepts name-value
%       arguments for program flow and visualisation options.
%
%   The function returns the updated table c with verdict and judged columns
%   populated. Save c to disk after the session to preserve your work.
%
% -------------------------------------------------------------------------
%   REQUIRED INPUTS
% -------------------------------------------------------------------------
%   c           - Detection table (captureHistoryTable). Must contain fields:
%                   t0, tEnd, detect_observer*, t0_observer*, tEnd_observer*,
%                   freq_observer*_1, freq_observer*_2, soundFolder_observer1
%                 On return, c will contain two additional fields:
%                   verdict  - 1 = true positive, 0 = false positive, NaN = unjudged
%                   judged   - logical flag, true if verdict has been assigned
%
%   sampleRate  - Target sample rate (Hz) for audio resampling and display.
%                 Audio is resampled to this rate before spectrogram computation.
%                 Example: 2000 (for low-frequency cetacean calls)
%
%   nfft        - FFT length (samples) for spectrogram computation.
%                 Larger values give finer frequency resolution but coarser
%                 time resolution. Must be a positive integer.
%                 Example: 512
%
%   noverlap    - Number of overlapping samples between adjacent FFT windows.
%                 Must be less than nfft. Higher overlap produces smoother
%                 spectrograms at the cost of computation time.
%                 Example: 256
%
% -------------------------------------------------------------------------
%   OPTIONAL POSITIONAL INPUTS
% -------------------------------------------------------------------------
%   pre         - Seconds of audio to display before the detection start (t0).
%                 Default: 30
%
%   post        - Seconds of audio to display after the detection end (tEnd).
%                 Default: 30
%
%   startIndex  - Row index in c at which to begin reviewing. Useful for
%                 resuming a session. Default: 1
%
%   increment   - Step size between detections. Use increment > 1 to review
%                 a subset (e.g., every 2nd detection). Default: 1
%
%   freqBand    - [low high] frequency limits (Hz) used to calculate the
%                 color scale percentiles. Only power within this band is
%                 used to set the display range, which helps focus contrast
%                 on the frequency band of interest.
%                 Default: [17 30]
%
%   suppressClicks - Name-value pair to suppress transient click noise:
%                   suppressClicks.threshold  Detection threshold (default: 3)
%                   suppressClicks.power      Suppression power (default: 1000)
%                   Pass nothing to disable click suppression (default).
%                   Example: judgeDetections(..., 'threshold', 3, 'power', 1000)
%
% -------------------------------------------------------------------------
%   OPTIONAL NAME-VALUE INPUTS — Program Flow
% -------------------------------------------------------------------------
%   showJudged      - If true, display detections that have already been
%                     judged (a banner shows the existing verdict). If false,
%                     already-judged detections are skipped silently.
%                     Default: false
%
%   overwriteJudged - If true, left/right clicks update the verdict even for
%                     detections already marked as judged. Has no effect
%                     unless showJudged is also true.
%                     Default: false
%
% -------------------------------------------------------------------------
%   OPTIONAL NAME-VALUE INPUTS — Visualisation
% -------------------------------------------------------------------------
%   windowType  - How to determine the audio window around each detection:
%                   'prePost'  Load pre seconds before t0 and post seconds
%                              after tEnd (default, adapts to call duration)
%                   'fixed'    Load a fixed windowDur-second window centred
%                              on the midpoint of the detection
%                   'tight'    Load only the detection itself (t0 to tEnd)
%                 Default: 'prePost'
%
%   windowDur   - Duration (s) of the fixed audio window. Only used when
%                 windowType is 'fixed'. Default: 60
%
%   cAdjust     - [floor ceiling] percentiles used for the colour scale,
%                 computed within freqBand. Values are clamped to [0, 100]
%                 and can be adjusted interactively with the arrow keys.
%                 Default: [5 99]
%
%   showObsBoxes - If true, overlay coloured bounding boxes for each
%                  observer's detection on the spectrogram. Each observer
%                  is assigned a distinct colour; the observer number is
%                  labelled at the top-left corner of each box.
%                  Default: true
%
%   cmPerPixel  - Width of each spectrogram time-bin in centimetres on
%                 screen. Controls how wide the figure window is drawn.
%                 Increase for more detail; decrease to fit more time on
%                 screen. Default: 0.05
%
%   figPos      - Starting position of the figure window in pixels,
%                 as [left bottom width height]. The width is overridden
%                 each frame by cmPerPixel, but left/bottom/height persist.
%                 Default: [100 400 679 400]
%
% -------------------------------------------------------------------------
%   CONTROLS
% -------------------------------------------------------------------------
%   VERDICT
%     Left click    Mark current detection as TRUE POSITIVE (verdict = 1)
%                   and advance to the next detection.
%     Right click   Mark current detection as FALSE POSITIVE (verdict = 0)
%                   and advance to the next detection.
%
%   INFORMATION
%     I             Print the capture history fields and values to console
%   AUDIO
%     P             Play the audio clip for the current detection window.
%   SCREENSHOT      
%     S             Save the adjudication window as a png named
%                   'adjudication_n.png where n is the row from the capture
%                   history table
%
%   DISPLAY — Arrow keys adjust the colour scale percentiles used for
%   contrast. The current floor/ceiling percentiles are shown in the title.
%
%     Left arrow    Lower the colour floor  -> darker minimum (more contrast)
%     Right arrow   Raise the colour floor  -> lighter minimum (less contrast)
%     Up arrow      Raise the colour ceiling -> lighter maximum
%     Down arrow    Lower the colour ceiling -> darker maximum
%
%   TIP: Press Left + Up to compress the colour range (more contrast).
%   Press Right + Down to expand it (less contrast).
%   Percentile limits are clamped to [0, 100].
%
% -------------------------------------------------------------------------
%   NOTES
% -------------------------------------------------------------------------
%   - The figure window width scales automatically with spectrogram length.
%   - If an error occurs during a detection, the function prints the stack
%     trace and exits, returning c as modified up to that point. Save c
%     immediately after a session ends.
%   - Audio is loaded using getAudioFromFiles, which must be
%     on the MATLAB path.
%
% -------------------------------------------------------------------------
%   EXAMPLES
% -------------------------------------------------------------------------
%   % Basic usage
%   c = judgeDetections(c, 2000, 512, 256);
%
%   % Custom pre/post window and frequency band
%   c = judgeDetections(c, 2000, 512, 256, 60, 60, 1, 1, [15 35]);
%
%   % Resume from row 50, reviewing every other detection
%   c = judgeDetections(c, 2000, 512, 256, 30, 30, 50, 2);
%
%   % Enable click suppression
%   c = judgeDetections(c, 2000, 512, 256, 30, 30, 1, 1, [17 30], ...
%       'threshold', 3, 'power', 1000);
%
%   % Review already-judged detections and allow overwriting
%   c = judgeDetections(c, 2000, 512, 256, 'showJudged', true, ...
%       'overwriteJudged', true);
%
%   % Use a fixed 90-second window and tighter initial colour scale
%   c = judgeDetections(c, 2000, 512, 256, 'windowType', 'fixed', ...
%       'windowDur', 90, 'cAdjust', [10 95]);
%
%   % Save results after session
%   save('myDetections.mat', 'c');
%
%   See also GINPUT, SPECTROGRAM, AUDIOPLAYER

% -------------------------------------------------------------------------
arguments
    c
    sampleRate   (1,1) {mustBeNumeric, mustBeNonempty}
    nfft         (1,1) {mustBeInteger, mustBePositive, mustBeNonempty}
    noverlap     (1,1) {mustBeInteger, mustBePositive, mustBeNonempty}
    pre          (1,1) {mustBeNumeric, mustBeNonempty} = 30 % seconds
    post         (1,1) {mustBeNumeric, mustBeNonempty} = 30 % seconds
    startIndex   (1,1) {mustBeInteger, mustBePositive} = 1
    increment    (1,1) {mustBeInteger, mustBePositive} = 1
    freqBand     (1,2) {mustBeNumeric}                 = [17 30] % Hz
    suppressClicks.threshold = []
    suppressClicks.power     = []
    % Program flow
    options.showJudged      (1,1) logical = false
    options.overwriteJudged (1,1) logical = false
    % Visualisation
    options.windowType  (1,1) string ...
        {mustBeMember(options.windowType,{'prePost','fixed','tight'})} = 'prePost'
    options.windowDur   (1,1) {mustBeNumeric, mustBePositive} = 60
    options.cAdjust     (1,2) {mustBeNumeric}                 = [5 99]
    options.showObsBoxes (1,1) logical                        = true
    options.cmPerPixel  (1,1) {mustBeNumeric, mustBePositive} = 0.05
    options.figPos      (1,4) {mustBeNumeric}                 = [100 400 679 400]
end

if (isempty(suppressClicks.threshold))
    suppressClicks = [];
end

% Unpack options into local variables for readability
showJudged      = options.showJudged;
overwriteJudged = options.overwriteJudged;
windowType      = options.windowType;
windowDur       = options.windowDur;
cAdjust         = options.cAdjust;
showObsBoxes    = options.showObsBoxes;
cmPerPixel      = options.cmPerPixel;
pos             = options.figPos;

%% Create verdict and judged columns unless they already exist
if ~any(strcmp('verdict',fieldnames(c)))
    c.verdict = nan(size(c.detect_observer1));
end
if ~any(strcmp('judged',fieldnames(c)))
    c.judged = false(size(c.verdict));
end

nDet = height(c);
wavInfo = wavFolderInfo(c.soundFolder_observer1{1});

figure('position',pos);
colormap(flipud(gray));
set(gcf,'units','centimeters');

i = startIndex;
while (i <= nDet)
    pos = get(gcf,'pos');
    msg = '';
    try
        if c.judged(i)==true && ~isnan(c.verdict(i))
            if showJudged
                msg = sprintf('Already judged as %g', c.verdict(i));
            else
                i = i+increment; % Uncomment to skip already judged
                continue
            end
        end % keep existing verdicts

        nObservers = sum(cellfun(@isempty,...
            strfind(c.Properties.VariableNames,'detect_observer') ...
            )==0);


        % calculate audio times
        callStart = min(c{i,append("t0_observer",string(1:nObservers))});
        callEnd = max(c{i,append("tEnd_observer",string(1:nObservers))});

        switch(windowType)
            case 'prePost'
                % Fixed pre and post
                startTime = callStart - pre/86400;
                endTime = callEnd + post/86400;
            case 'fixed'
                startTime = mean([callStart callEnd])-0.5*windowDur/86400;
                endTime = mean([callStart callEnd])+0.5*windowDur/86400;
            otherwise
                startTime = callStart;
                endTime = callEnd;
        end
        wav =  getAudioFromFiles(wavInfo,startTime,endTime,newRate=sampleRate);

        % Remove clicks
        if ~isempty(suppressClicks)
            threshold = 3;
            power = 1000;
            if isfield(suppressClicks,'threshold')
                threshold = suppressClicks.threshold;
            end
            if isfield(suppressClicks,'power')
                power = suppressClicks.power;
            end
            wav = removeClicks(wav,threshold,power);
        end



        [s,f,t,p] = spectrogram(wav,nfft,noverlap,[],sampleRate,'yAxis');

        pdB = 20*log10(p);
        fIx = double(f > freqBand(1) & f < freqBand(2));
        nanIx = find(fIx==0);
        fIx(nanIx)=nan(size(nanIx));
        pIx = repmat(fIx,1,length(t));
        %         h = pcolor(t,f,pdB.*pIx);
        h = pcolor(t,f,pdB);

        set(h,'lineStyle','none');
        %         view([0,90]);
        %         clim(1) = min(pdB(:).*pIx(:));
        clim(1) = prctile(pdB(:).*pIx(:),cAdjust(1));
        %         clim(2) = max(pdB(:).*pIx(:));
        clim(2) = prctile(pdB(:).*pIx(:),cAdjust(2));

        set(gca,'clim',clim);

        %         caxis([100,200]-200);
        %         ylim([0,0.5]);
        horzline(freqBand,'c--');
        vertline([(callStart-startTime)*86400 (callEnd-startTime)*86400],'c--')


        if showObsBoxes
            l = c{i,append("t0_observer",string(1:nObservers))}-startTime;
            r = c{i,append("tEnd_observer",string(1:nObservers))}-startTime;
            try % Format as of 2026-03-18
                top = c{i,append("freq_observer",string(1:nObservers),'_1')};
                bot = c{i,append("freq_observer",string(1:nObservers),'_2')};
            catch % Format from 2025-10-31
                top = c{i,append("freq_1_observer",string(1:nObservers))};
                bot = c{i,append("freq_2_observer",string(1:nObservers))};
            end
            ho = ishold;
            hold on;
            obsColors = lines(nObservers);
            for j = 1:nObservers
                plotBox([l(j),r(j)]*86400,[bot(j),top(j)],obsColors(j,:),'linewidth',1.5);
                text(l(j)*86400,top(j),num2str(j),'color',obsColors(j,:))
            end
            if ho==0; hold off; end
        end

        ti = {sprintf('#% 4g/%g floor/ceil percentile=[%g %g]',...
            1+(i-1)/increment,floor(nDet/increment), cAdjust), ...
            sprintf('(Row % 4g: %s)',i, datestr(c.t0(i)))};

        title(ti);
        text(mean(xlim),max(ylim),msg,'VerticalAlignment','top',...
            'FontSize',14,'FontWeight','bold','HorizontalAlignment','center','BackgroundColor','w');

        pos(3) = size(s,2)*cmPerPixel;
        set(gcf,'position',pos);
        drawnow;

        [~,~,button] = ginput(1);
        switch button
            case 28 % Left arrow shift color floor down (darker)
                cAdjust = cAdjust+[-1 0];
            case 29
                cAdjust = cAdjust+[1 0]; % Down arrow color floor down (ligher)
            case 30 % Up arrow shift color ceiling up (lighter)
                cAdjust = cAdjust+[0 1];
            case 31
                cAdjust = cAdjust+[0 -1]; % Down arrow to shift cieling down (darker)
            case 1 % Left click for true positive;
                % Only write if new verdicts (unless overwrite is true)
                if c.judged(i)~=true || overwriteJudged
                    c.verdict(i) = 1;
                    c.judged(i)=true;
                end
                i = i + increment;
            case 3 % Right click for false positive
                % Only write new verdicts or if overwriteJudged==true)
                if c.judged(i)~=true || overwriteJudged
                    c.verdict(i) = 0;
                    c.judged(i)=true;
                end
                i = i + increment;
            case {2,45} % -/Minus-key/Middle click for 'other'
                % Only write if new verdicts (unless overwrite is true)
                if c.judged(i)~=true || overwriteJudged
                    c.verdict(i) = -1;
                    c.judged(i)=true;
                end
                i = i + increment;
            case 105 % i key - print capture history row in console
                disp(table2struct(c(i,:)))
            case 112 % p key - Play visible audio
                player = audioplayer(wav,sampleRate);
                play(player);
            case 115 % s key - save screenshot
                saveAsPng(sprintf('adjudication_%03g',i))
                
        end

        % Ensure percentile stay within 0 and 100;
        cAdjust(1) = max(cAdjust(1),0);
        cAdjust(2) = min(cAdjust(2),100);

        disp(ti);
    catch
        l = lasterror;
        disp(l);
        for j = 1:length(l.stack); disp(l.stack(j)); end
        fprintf('Exiting on index %g\n',i);
        return
    end
end
