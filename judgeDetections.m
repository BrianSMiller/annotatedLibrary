function c = judgeDetections(c,sampleRate,nfft,noverlap,pre,post, ...
    startIndex,increment,freqBand,options)

% JUDGEDETECTIONS  Manually review and label candidate detections
%
%   c = judgeDetections(c, sampleRate, nfft, noverlap) launches an
%   interactive spectrogram viewer for manually labeling detections in a
%   capture history table as true positives or false positives.
%
%   c = judgeDetections(c, sampleRate, nfft, noverlap, pre, post, ...
%       startIndex, increment, freqBand) uses additional positional options
%       described below.
%
%   c = judgeDetections(___, Name, Value) additionally accepts name-value
%       arguments for click suppression, program flow, and visualisation.
%
%   The function returns the updated table c with verdict and judged columns
%   populated. Save c to disk after the session to preserve your work.
%
% -------------------------------------------------------------------------
%   REQUIRED INPUTS
% -------------------------------------------------------------------------
%   c           - Detection table (captureHistoryTable). Must contain fields:
%                   t0, tEnd, detect_observer*, t0_observer*, tEnd_observer*,
%                   soundFolder_observer1 (or soundFolder; see NOTES).
%                 For showObsBoxes, also needs per-observer frequency bounds
%                 under one of the naming conventions in NOTES.
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
% -------------------------------------------------------------------------
%   OPTIONAL NAME-VALUE INPUTS — Click Suppression
% -------------------------------------------------------------------------
%   suppressClickThreshold - Detection threshold for transient click
%                            suppression. Higher values suppress only
%                            louder clicks. Set to [] to disable.
%                            Default: [] (suppression disabled)
%
%   suppressClickPower     - Suppression power applied to detected clicks.
%                            Higher values more aggressively attenuate.
%                            Only used if suppressClickThreshold is set.
%                            Default: [] (suppression disabled)
%
%                   Example: judgeDetections(c, 2000, 512, 256, ...
%                       'suppressClickThreshold', 3, 'suppressClickPower', 1000)
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
%   soundFolder - Explicit path to the folder of wav files, overriding
%                 whatever the capture history table would otherwise
%                 resolve to. Use this when you want to bypass table-based
%                 resolution entirely. See NOTES for the full resolution
%                 order.
%                 Default: [] (resolve from c; see NOTES)
%
%   showObsBoxes - If true, overlay coloured bounding boxes for each
%                  observer's detection on the spectrogram. Each observer
%                  is assigned a distinct colour; the observer number is
%                  labelled at the top-left corner of each box.
%                  Default: true
%
%   showRefinedBounds - If true, overlay a magenta box using c's
%                  top-level t0/tEnd/fLow/fHigh -- whatever bounds are
%                  CURRENTLY primary for each detection (refined in
%                  place if a bounds-refining SNR method like kalmanRidge
%                  genuinely tracked a call). If c has a boundsRefined
%                  column and it's false for a row (dropout -- t0/tEnd/
%                  freq pass through the ORIGINAL bounds unchanged per
%                  bsnr's contract), the box is skipped and an explicit
%                  "DROPOUT" flag is shown instead, since drawing a
%                  "refined" box there would land exactly on top of the
%                  original-bounds box and look like a rendering bug, not
%                  the real underlying signal. If boundsRefined isn't
%                  present at all, always draws the box (no way to tell
%                  refined from pass-through in that case).
%                  Requires top-level fLow/fHigh columns; if absent,
%                  disabled with a one-time console note rather than
%                  erroring per detection.
%                  Default: true
%
%   showOriginalBounds - If true, overlay a white dashed box using
%                  c.originalT0/originalTEnd/originalFLow/originalFHigh
%                  -- the detector's bounds BEFORE any refinement. Always
%                  drawn when available, regardless of dropout status.
%                  Requires all four columns; if absent (e.g. c never
%                  went through a bounds-refining SNR step), disabled
%                  with a one-time console note.
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
%     Left click          Mark current detection as TRUE POSITIVE (verdict = 1)
%                         and advance to the next detection.
%     Right click         Mark current detection as FALSE POSITIVE (verdict = 0)
%                         and advance to the next detection.
%     Middle click or -   Mark current detection as OTHER (verdict = -1)
%                         and advance to the next detection.
%
%   AUDIO
%     P                   Play the audio clip for the current detection window.
%
%   INFORMATION
%     I                   Print the capture history fields and values for the
%                         current detection to the console.
%
%   SCREENSHOT
%     S                   Save the current figure as 'adjudication_NNN.png'
%                         where NNN is the row index in c.
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
%   - soundFolder is resolved ONCE per session, before the review loop,
%     and reused for every row. c is assumed to be a single-deployment
%     table. Resolution order:
%       1. options.soundFolder, if given.
%       2. c.soundFolder, if the table has a top-level column (e.g. from
%          resolveObserverColumns).
%       3. First non-empty soundFolder_observerK across all rows and all
%          observers (legacy fallback for tables predating
%          resolveObserverColumns).
%     Errors immediately, before opening a figure, if none resolve, or if
%     c.soundFolder (step 2) holds more than one distinct value -- a
%     mixed-deployment table needs to be split before review, since only
%     one soundFolder is ever loaded.
%   - Per-observer frequency bounds for showObsBoxes are resolved once,
%     before the review loop, checking (in order): fHigh_observerK /
%     fLow_observerK (current matchbox convention), freq_observerK_1/_2
%     (2026-03-18), freq_1_observerK/freq_2_observerK (2025-10-31). Errors
%     immediately if none match.
%   - showRefinedBounds/showOriginalBounds availability is checked once,
%     before the review loop, not per detection. If required columns are
%     missing, the option is silently disabled with a one-time console
%     note rather than erroring on the first detection. showRefinedBounds
%     needs top-level t0/tEnd/fLow/fHigh; showOriginalBounds needs
%     originalT0/originalTEnd/originalFLow/originalFHigh, which only
%     exist if c went through a bounds-refining SNR step (e.g. qcAttachSnr
%     with snrType='kalmanRidge' or 'momentumRidge'). If c also has a
%     boundsRefined column, showRefinedBounds uses it PER ROW: a false
%     value (dropout) shows an orange "DROPOUT" flag instead of a box,
%     rather than drawing a refined box that would coincide exactly with
%     the original-bounds box and look broken. Without a boundsRefined
%     column, the refined box is always drawn (no way to distinguish
%     genuine refinement from pass-through in that case).
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
%   c = judgeDetections(c, 2000, 512, 256, 'suppressClickThreshold', 3, ...
%       'suppressClickPower', 1000);
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
    % Click suppression
    options.suppressClickThreshold = []
    options.suppressClickPower     = []
    % Program flow
    options.showJudged      (1,1) logical = false
    options.overwriteJudged (1,1) logical = false
    % Visualisation
    options.windowType  (1,1) string ...
        {mustBeMember(options.windowType,{'prePost','fixed','tight'})} = 'prePost'
    options.windowDur   (1,1) {mustBeNumeric, mustBePositive} = 60
    options.cAdjust     (1,2) {mustBeNumeric}                 = [5 99]
    options.soundFolder  {mustBeText}                         = ''
    options.showObsBoxes (1,1) logical                        = true
    options.showRefinedBounds  (1,1) logical                  = true
    options.showOriginalBounds (1,1) logical                  = true
    options.cmPerPixel  (1,1) {mustBeNumeric, mustBePositive} = 0.05
    options.figPos      (1,4) {mustBeNumeric}                 = [100 400 679 400]
end

% Unpack options into local variables for readability
showJudged      = options.showJudged;
overwriteJudged = options.overwriteJudged;
windowType      = options.windowType;
windowDur       = options.windowDur;
cAdjust         = options.cAdjust;
showObsBoxes    = options.showObsBoxes;
showRefinedBounds  = options.showRefinedBounds;
showOriginalBounds = options.showOriginalBounds;
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
soundFolder = resolveSoundFolder(c, options.soundFolder);
wavInfo = wavFolderInfo(soundFolder);

% Exact match on detect_observerN only -- a naive substring search
% (strfind) also catches double-suffixed columns like
% detect_observer1_observer1, which appear when matchbox is used to
% match an already-matched CHT against a further detection table (e.g.
% production output vs ATBFL): matchbox blindly re-suffixes every column
% of each input, including columns that already carry their own
% detect_observerK suffix from a prior match. That inflates nObservers
% (5 instead of 2 in the case that surfaced this) and breaks every
% downstream lookup that assumes detect_observer1..detect_observerN with
% no gaps. Anchored regex ('^...$') only matches the genuine, currently-
% suffixed columns, regardless of how many prior matching passes a
% column's name has already been through.
nObservers = sum(~cellfun(@isempty, ...
    regexp(c.Properties.VariableNames, '^detect_observer\d+$', 'once') ...
    ));
[topFreqCols, botFreqCols] = resolveObserverFreqColumns(c, nObservers);

if showRefinedBounds && ~all(ismember({'t0','tEnd','fLow','fHigh'}, c.Properties.VariableNames))
    fprintf(['judgeDetections: showRefinedBounds requested but c is missing top-level ' ...
        't0/tEnd/fLow/fHigh -- disabling.\n']);
    showRefinedBounds = false;
end

if showOriginalBounds && ~all(ismember( ...
        {'originalT0','originalTEnd','originalFLow','originalFHigh'}, c.Properties.VariableNames))
    fprintf(['judgeDetections: showOriginalBounds requested but c is missing ' ...
        'originalT0/originalTEnd/originalFLow/originalFHigh (only present if c went ' ...
        'through a bounds-refining SNR step, e.g. qcAttachSnr with kalmanRidge) -- disabling.\n']);
    showOriginalBounds = false;
end

hasBoundsRefinedCol = any(strcmp('boundsRefined', c.Properties.VariableNames));

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

        % nObservers resolved once, before the loop (see above)


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
        if ~isempty(options.suppressClickThreshold)
            threshold = options.suppressClickThreshold;
            power = 1000;
            if ~isempty(options.suppressClickPower)
                power = options.suppressClickPower;
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
        yline(freqBand,'c--');
        xline([(callStart-startTime)*86400 (callEnd-startTime)*86400],'c--')


        if showObsBoxes
            l = c{i,append("t0_observer",string(1:nObservers))}-startTime;
            r = c{i,append("tEnd_observer",string(1:nObservers))}-startTime;
            top = c{i, topFreqCols};
            bot = c{i, botFreqCols};
            ho = ishold;
            hold on;
            obsColors = lines(nObservers);
            for j = 1:nObservers
                plotBox([l(j),r(j)]*86400,[bot(j),top(j)],obsColors(j,:),'linewidth',1.5);
                text(l(j)*86400,top(j),num2str(j),'color',obsColors(j,:))
            end
            if ho==0; hold off; end
        end

        if showRefinedBounds
            isDropout = hasBoundsRefinedCol && ~c.boundsRefined(i);
            ho = ishold;
            hold on;
            if isDropout
                % On dropout, t0/tEnd/freq pass through the ORIGINAL bounds
                % unchanged (bsnr's documented contract) -- drawing a
                % "refined" box here would land exactly on top of the
                % original-bounds box, indistinguishable from it. Flag it
                % explicitly instead of drawing a coincident, misleading box.
                text(mean(xlim),min(ylim),'kalmanRidge: DROPOUT (no refinement -- bounds = original)', ...
                    'color',[1 0.6 0],'FontSize',9,'FontWeight','bold', ...
                    'VerticalAlignment','bottom','HorizontalAlignment','center','BackgroundColor','k');
            else
                l = (c.t0(i) - startTime)*86400;
                r = (c.tEnd(i) - startTime)*86400;
                plotBox([l,r],[c.fLow(i),c.fHigh(i)],'m','linewidth',2);
                text(l,c.fHigh(i),'refined','color','m','FontSize',8);
            end
            if ho==0; hold off; end
        end

        if showOriginalBounds
            ho = ishold;
            hold on;
            l = (c.originalT0(i) - startTime)*86400;
            r = (c.originalTEnd(i) - startTime)*86400;
            plotBox([l,r],[c.originalFLow(i),c.originalFHigh(i)],'w','linewidth',1,'linestyle','--');
            text(l,c.originalFHigh(i),'orig','color','w','FontSize',8);
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
end

function soundFolder = resolveSoundFolder(c, override)
% resolveSoundFolder  Priority-ordered soundFolder resolution for
% judgeDetections. See judgeDetections' NOTES section for the resolution
% order this implements. Errors immediately (before any figure is opened)
% rather than letting wavFolderInfo fail on an empty/ambiguous path.

% 1. Explicit override always wins.
if ~isempty(override)
    soundFolder = char(override);
    return
end

% 2. Resolved top-level column (e.g. from resolveObserverColumns).
if any(strcmp('soundFolder', c.Properties.VariableNames))
    vals = c.soundFolder;
    nonEmpty = vals(~cellfun(@isempty, vals));
    if isempty(nonEmpty)
        error('judgeDetections:emptySoundFolderColumn', ...
            'c.soundFolder exists but every entry is empty.');
    end
    uVals = unique(nonEmpty);
    if numel(uVals) > 1
        error('judgeDetections:mixedSoundFolder', ...
            ['c.soundFolder has %d distinct values. judgeDetections loads audio ' ...
             'from ONE folder for the whole session -- split this table by ' ...
             'deployment first, or pass ''soundFolder'' explicitly.'], numel(uVals));
    end
    soundFolder = uVals{1};
    return
end

% 3. Legacy fallback: first non-empty value across all soundFolder_observerK
%    columns and all rows (not just observer1, row1).
obsCols = c.Properties.VariableNames(startsWith(c.Properties.VariableNames, 'soundFolder_observer'));
for k = 1:numel(obsCols)
    vals = c.(obsCols{k});
    nonEmpty = vals(~cellfun(@isempty, vals));
    if ~isempty(nonEmpty)
        soundFolder = nonEmpty{1};
        return
    end
end

% 4. Nothing resolved -- fail loud, before opening a figure.
error('judgeDetections:noSoundFolder', ...
    ['Could not resolve a soundFolder from c. Pass ''soundFolder'' explicitly, ' ...
     'or run resolveObserverColumns(c, {''soundFolder''}) first.']);
end

function [topCols, botCols] = resolveObserverFreqColumns(c, nObservers)
% resolveObserverFreqColumns  Find the per-observer high/low frequency
% column names actually present in c, across the naming conventions
% matchbox and judgeDetections have used over time. Checked ONCE before
% the review loop starts, rather than guessed per-detection via try/catch.

vn = c.Properties.VariableNames;

% Current matchbox convention (fLow/fHigh -- see multiCaptureHistory*.m)
topA = append("fHigh_observer", string(1:nObservers));
botA = append("fLow_observer",  string(1:nObservers));
if all(ismember(topA, vn)) && all(ismember(botA, vn))
    topCols = topA; botCols = botA; return
end

% Format as of 2026-03-18
topB = append("freq_observer", string(1:nObservers), '_1');
botB = append("freq_observer", string(1:nObservers), '_2');
if all(ismember(topB, vn)) && all(ismember(botB, vn))
    topCols = topB; botCols = botB; return
end

% Format from 2025-10-31
topC = append("freq_1_observer", string(1:nObservers));
botC = append("freq_2_observer", string(1:nObservers));
if all(ismember(topC, vn)) && all(ismember(botC, vn))
    topCols = topC; botCols = botC; return
end

error('judgeDetections:noFreqColumns', ...
    ['Could not find per-observer high/low frequency columns under any known ' ...
     'naming convention (fHigh_observerK/fLow_observerK, freq_observerK_1/_2, ' ...
     'or freq_1_observerK/freq_2_observerK). Check c.Properties.VariableNames.']);
end
