function detTable = qcValidateDetections(detTable, siteMeta, varargin)
% qcValidateDetections  Simple heuristic sanity checks on a detection
% table, flagging physically or logistically implausible detections
% before they reach SNR estimation or visual QC.
%
% Adds two columns (does not remove any rows -- downstream steps decide
% what to do with qcValid==false):
%   qcValid - true if the detection passes all enabled checks
%   qcFlag  - string naming the first failing check ('' if valid)
%
% Checks (each toggleable via checkXxx, true by default):
%   checkIncreasingFreq   - fLow < fHigh, neither NaN
%   checkMinFrequency     - fLow/fHigh finite and fLow > minFreq_Hz
%                           (default 5 Hz -- comfortably below real whale
%                           calls, but catches near-zero/normalized
%                           frequency values that never got converted to Hz)
%   checkMinBandwidth     - fHigh-fLow > minBandwidth_Hz (default 2 Hz) --
%                           near-zero bandwidth is implausible for any
%                           real frequency-swept call, independent of the
%                           specific species/call-type's true bandwidth
%   checkNyquist          - fHigh < actual audio sample rate / 2, read
%                           from wavFolderInfo(siteMeta.wavFolder) rather
%                           than siteMeta.calib.sampleRate -- the latter
%                           may describe a different stage of the
%                           pipeline (e.g. original recorder config vs. a
%                           decimated product actually on disk)
%   checkDuration         - tEnd-t0 within [minDuration_s maxDuration_s]
%                           (default 0.5-60 s)
%   checkDeploymentBounds - t0/tEnd within siteMeta.calib.startDate/endDate.
%                           Violations within boundsClampTolerance_s are
%                           clamped to the boundary and kept (flagged
%                           "ClampedToDeploymentWindow" but qcValid stays
%                           true); larger violations are flagged invalid
%                           ("OutsideDeploymentWindow"). A few ms/samples
%                           over is usually an off-by-one-frame or
%                           rounding issue; way outside the window is a
%                           much more fundamental problem worth treating
%                           the whole batch with caution over.
%   checkWavBounds        - same clamp-vs-invalid logic as
%                           checkDeploymentBounds, but against a SINGLE
%                           file's own [startDate, startDate+duration]
%                           window using per-file data from
%                           wavFolderInfo(siteMeta.wavFolder), rather than
%                           an overall folder min/max. This matters
%                           because sample rate can drift and different
%                           software (Raven, PAMGuard) disagrees on
%                           whether within-file offsets are counted from
%                           file start or calendar time -- per-file
%                           containment catches both a detection sitting
%                           in a genuine gap and one that's drifted across
%                           a file boundary, neither of which an overall
%                           bounds check would catch.
%   checkDuplicates       - exact repeated (t0,tEnd,fLow,fHigh) rows
%                           (possible double-counted detections)
%
%   detTable = qcValidateDetections(detTable, siteMeta)
%   detTable = qcValidateDetections(detTable, siteMeta, ...
%                   'minDuration_s', 0.25, 'minFreq_Hz', 8, ...
%                   'checkDuplicates', false)
%
% Inputs
%   detTable  - detection table with t0, tEnd (datenum), fLow/fHigh (or
%               freq), soundFolder
%   siteMeta  - the site's metadata struct directly, e.g.
%               metaDataKerguelen2015_250Hz() -- .wavFolder and
%               .startDate/.endDate (datenum) used directly, no nested
%               .calib wrapper. Nyquist reads the actual sample rate from
%               wavFolderInfo(siteMeta.wavFolder) rather than any
%               sampleRate field here (see check 3 below for why). Any
%               check whose required field is missing is skipped with a
%               warning rather than erroring.
%
% B. Miller, AAD, 2026

p = inputParser;
addParameter(p,'checkWavBounds',true);
addParameter(p,'checkDeploymentBounds',true);
addParameter(p,'checkNyquist',true);
addParameter(p,'checkMinFrequency',true);
addParameter(p,'checkIncreasingFreq',true);
addParameter(p,'checkMinBandwidth',true);
addParameter(p,'minBandwidth_Hz',2);
addParameter(p,'checkDuration',true);
addParameter(p,'checkDuplicates',true);
addParameter(p,'minDuration_s',0.128);
addParameter(p,'maxDuration_s',60);
addParameter(p,'minFreq_Hz',5);
addParameter(p,'boundsClampTolerance_s',1);
parse(p,varargin{:});
opt = p.Results;

nDet = height(detTable);

if ~any(strcmp('fLow',detTable.Properties.VariableNames)) && ...
        any(strcmp('freq',detTable.Properties.VariableNames))
    detTable.fLow  = detTable.freq(:,1);
    detTable.fHigh = detTable.freq(:,2);
end

qcValid = true(nDet,1);
qcFlag  = repmat("",nDet,1);

%% Fetch real audio properties once, used by both Nyquist and wavBounds
% checks. Deliberately NOT using siteMeta.calib.sampleRate here -- for
% Kerguelen2015 (and likely others) that field reflects the ORIGINAL
% recorder configuration (12000 Hz), not the actual audio sitting in
% wavFolder (a decimated 250 Hz product). Reading the real per-file rate
% from wavFolderInfo avoids trusting a metadata field that may describe a
% different stage of the pipeline than the audio actually being analyzed.
wi = [];
if (opt.checkWavBounds || opt.checkNyquist) && ...
        isfield(siteMeta,'wavFolder') && ~isempty(siteMeta.wavFolder)
    try
        % NOTE: passing [] not '' for timeStampFormat -- the pasted
        % wavFolderInfo docstring's own silent-mode example uses '', but
        % that contradicts the established rule that '' hangs. Trusting
        % the hard-won rule over the docstring example here; worth
        % confirming which is actually still true.
        wi = wavFolderInfo(siteMeta.wavFolder, [], false, false);
    catch ME
        warning('qcValidateDetections:wavInfoFailed', ...
            'Could not read wavFolderInfo(%s): %s -- Nyquist/wavBounds checks skipped.', ...
            siteMeta.wavFolder, ME.message);
    end
end

%% 1. Increasing frequency range
if opt.checkIncreasingFreq
    fail = isnan(detTable.fLow) | isnan(detTable.fHigh) | ...
           (detTable.fHigh <= detTable.fLow);
    qcFlag(fail & qcValid) = "NonIncreasingFreq";
    qcValid(fail) = false;
end

%% 2. Minimum plausible frequency (also catches non-finite values)
if opt.checkMinFrequency
    fail = ~isfinite(detTable.fLow) | ~isfinite(detTable.fHigh) | ...
           detTable.fLow <= opt.minFreq_Hz;
    qcFlag(fail & qcValid) = "ImplausibleFrequency";
    qcValid(fail) = false;
end

%% 2b. Minimum plausible bandwidth
% A near-zero fHigh-fLow bandwidth is implausible for any real
% frequency-swept call (D-call, fin 40Hz downsweep, etc.) regardless of
% exact species/call-type expectations -- more likely a narrowband
% transient or noise blip misclassified as a detection. Conservative
% default (2 Hz) deliberately doesn't assume a specific call's true
% bandwidth; tighten if you know the target call's expected bandwidth.
if opt.checkMinBandwidth
    bw = detTable.fHigh - detTable.fLow;
    fail = qcValid & (bw < opt.minBandwidth_Hz);
    qcFlag(fail & qcValid) = "ImplausibleBandwidth";
    qcValid(fail) = false;
end

%% 3. Nyquist -- uses real sample rate from wi, not siteMeta.calib.sampleRate
if opt.checkNyquist
    if ~isempty(wi)
        rates = unique([wi.sampleRate]);
        if numel(rates) > 1
            warning('qcValidateDetections:mixedSampleRates', ...
                ['Audio in %s reports %d different sample rates (%s Hz) -- ' ...
                 'using the minimum for the Nyquist check (most conservative). ' ...
                 'Worth checking whether the recorder configuration changed ' ...
                 'mid-deployment.'], siteMeta.wavFolder, numel(rates), ...
                strtrim(sprintf('%g ', rates)));
        end
        nyquistHz = min(rates)/2;
        fail = detTable.fHigh >= nyquistHz;
        qcFlag(fail & qcValid) = "ExceedsNyquist";
        qcValid(fail) = false;
    else
        warning('qcValidateDetections:noSampleRate', ...
            'Could not determine actual audio sample rate -- skipping Nyquist check.');
    end
end

%% 4. Duration sanity
if opt.checkDuration
    durSec = (detTable.tEnd - detTable.t0) * 86400;
    fail = durSec < opt.minDuration_s | durSec > opt.maxDuration_s;
    qcFlag(fail & qcValid) = "ImplausibleDuration";
    qcValid(fail) = false;
end

%% 5. Deployment bounds
if opt.checkDeploymentBounds
    if all(isfield(siteMeta,{'startDate','endDate'}))
        beforeStartDays = max(siteMeta.startDate - detTable.t0, 0);
        afterEndDays    = max(detTable.tEnd - siteMeta.endDate, 0);
        violationSec    = max(beforeStartDays, afterEndDays) * 86400;

        [clampable, severe] = classifyBoundsViolation(violationSec, ...
            opt.boundsClampTolerance_s, 'Deployment-window');

        if any(clampable)
            clampStart = clampable & beforeStartDays > 0;
            clampEnd   = clampable & afterEndDays > 0;
            detTable.t0(clampStart)  = siteMeta.startDate;
            detTable.tEnd(clampEnd)  = siteMeta.endDate;
            qcFlag(clampable & qcValid) = "ClampedToDeploymentWindow";
        end

        qcFlag(severe & qcValid) = "OutsideDeploymentWindow";
        qcValid(severe) = false;
    else
        warning('qcValidateDetections:noDeploymentWindow', ...
            'siteMeta.startDate/endDate not found -- skipping deployment-bounds check.');
    end
end

%% 6. Per-file wav bounds (not just overall folder min/max)
if opt.checkWavBounds
    if ~isempty(wi)
        try
            % Root cause of the previously-seen "different number of
            % elements" error, now confirmed: [wi.startDate] (and hence
            % `ends`) comes out as a ROW vector (struct-array bracket
            % concatenation always produces a row), while detTable.t0/
            % tEnd are COLUMNs. Indexing a vector preserves ITS OWN
            % orientation regardless of the index array's shape, so
            % ends(fileIx) stayed a row even with fileIx as a column --
            % column minus row then silently broadcasts into a matrix
            % instead of subtracting elementwise, which fails on
            % assignment back into a column. Forcing column orientation
            % throughout fixes this.
            startsCol = [wi.startDate];
            startsCol = startsCol(:);
            dursCol = [wi.duration];
            dursCol = dursCol(:);
            [startsCol, sortIx] = sort(startsCol);
            dursCol = dursCol(sortIx);
            endsCol = startsCol + dursCol/86400;

            t0col   = detTable.t0(:);
            tEndCol = detTable.tEnd(:);

            fileIx = discretize(t0col, [startsCol; Inf]);

            % t0 slightly before the very first file also gets a chance
            % to clamp rather than being an automatic "no file matched"
            % case -- same off-by-a-few-samples reasoning as the
            % deployment-window check, just at the very start of the
            % recording rather than the deployment overall.
            beforeFirstDays = max(startsCol(1) - t0col, 0);
            nearFirst = isnan(fileIx) & (beforeFirstDays*86400 <= opt.boundsClampTolerance_s);
            fileIx(nearFirst) = 1;

            validFile = ~isnan(fileIx);

            startViolDays = zeros(nDet,1);
            startViolDays(nearFirst) = beforeFirstDays(nearFirst);

            endViolDays = zeros(nDet,1);
            endViolDays(validFile) = max(tEndCol(validFile) - endsCol(fileIx(validFile)), 0);

            violationDays = max(startViolDays, endViolDays);
            violationDays(~validFile) = Inf; % genuinely unmatched to any file

            violationSec = violationDays * 86400;

            [clampable, severe] = classifyBoundsViolation(violationSec, ...
                opt.boundsClampTolerance_s, 'Audio-bounds');

            if any(clampable)
                clampStart = clampable & startViolDays > 0;
                clampEnd   = clampable & endViolDays > 0;
                detTable.t0(clampStart)  = startsCol(1);
                detTable.tEnd(clampEnd)  = endsCol(fileIx(clampEnd));
                qcFlag(clampable & qcValid) = "ClampedToAudioBounds";
            end

            qcFlag(severe & qcValid) = "OutsideAudioBounds";
            qcValid(severe) = false;
        catch ME
            warning('qcValidateDetections:wavBoundsFailed', ...
                'Could not process per-file audio bounds -- skipping wav-bounds check.\n%s', ...
                getReport(ME, 'extended', 'hyperlinks', 'off'));
        end
    else
        warning('qcValidateDetections:noWavFolder', ...
            'Could not determine audio bounds -- skipping wav-bounds check.');
    end
end

%% 7. Duplicate detections
if opt.checkDuplicates
    key = [detTable.t0, detTable.tEnd, detTable.fLow, detTable.fHigh];
    [~, firstIx] = unique(key,'rows','stable');
    isDup = true(nDet,1);
    isDup(firstIx) = false;
    qcFlag(isDup & qcValid) = "DuplicateDetection";
    qcValid(isDup) = false;
end

detTable.qcValid = qcValid;
detTable.qcFlag  = qcFlag;

nBad = sum(~qcValid);
if nBad > 0
    fprintf('qcValidateDetections: %d of %d detections flagged:\n', nBad, nDet);
    disp(groupsummary(detTable(~qcValid,:),'qcFlag'));
end

end

function [clampable, severe] = classifyBoundsViolation(violationSec, toleranceSec, label)
% Splits a vector of time-bounds violations (seconds, 0 = no violation,
% Inf = no matching reference at all) into clampable (small, within
% tolerance -- likely a rounding/off-by-one-frame issue) vs severe
% (beyond tolerance -- a real data-quality problem). Prints magnitude
% stats so the person can judge severity at a glance rather than just
% getting a bare count.
fail = violationSec > 0;
clampable = fail & isfinite(violationSec) & (violationSec <= toleranceSec);
severe    = fail & ~clampable;

if ~any(fail)
    return;
end

finiteViol = violationSec(fail & isfinite(violationSec));
nInf = sum(isinf(violationSec) & fail);
extra = '';
if nInf > 0
    extra = sprintf(', plus %d detection(s) with no matching reference at all', nInf);
end

if ~isempty(finiteViol)
    fprintf(['qcValidateDetections: %s violations for %d detection(s) -- ' ...
        'min %.3f s, median %.3f s, max %.3f s%s ' ...
        '(clamp tolerance: %.3f s; %d clamped, %d flagged invalid)\n'], ...
        label, sum(fail), min(finiteViol), median(finiteViol), max(finiteViol), extra, ...
        toleranceSec, sum(clampable), sum(severe));
else
    fprintf(['qcValidateDetections: %s violations for %d detection(s)%s ' ...
        '(clamp tolerance: %.3f s; %d clamped, %d flagged invalid)\n'], ...
        label, sum(fail), extra, toleranceSec, sum(clampable), sum(severe));
end
end
