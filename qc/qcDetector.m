function d = qcDetector(cfg, useCalibration)
% qcDetector  Run the full QC pipeline for a single detector config,
% single site: load+convert -> validate -> attach SNR (valid rows only)
% -> diagnostic plots -> full-set visual paging.
%
% Deliberately NOT generalized to multi-site or filename-based site
% inference -- there's one real site (Kerguelen2015) today, detector
% output/format changes constantly, and building machinery for
% hypothetical future sites isn't worth it yet. When a detector actually
% covers multiple sites, extend this then, with real requirements in hand.
%
%   d = qcDetector(cfg, useCalibration)
%
% Inputs
%   cfg  - struct with fields:
%            name          - detector label, e.g. 'Coni_LFDCS_Dcall'
%            folder        - folder containing the detection file
%            filePattern   - exact filename to load
%            converter     - function handle, e.g. @lfdcsTableToDetection
%            siteCode      - site code, e.g. 'Kerguelen2015' (used for
%                            the 'site' column and diagnostics only)
%            siteMeta      - metadata struct for this exact site/product,
%                            e.g. metaDataKerguelen2015_250Hz() -- caller
%                            resolves this directly so there's no
%                            ambiguity about which metaData function is
%                            the right one for the audio actually in use
%            snrType       - passed to qcAttachSnr
%            nfft          - passed to qcAttachSnr
%            spectroParams - passed to qcPageDetections
%   useCalibration - if true, attach calibrated SNR using cfg.siteMeta;
%                    if false, uncalibrated RMS
%
% Output
%   d  - full detection table (valid + flagged rows), post-SNR-attach,
%        with qcValid/qcFlag columns from qcValidateDetections, and (as
%        of 2026-07-24) t0/tEnd/freq possibly REFINED in place by bsnr
%        (momentumRidge/kalmanRidge) -- see originalT0/originalTEnd/
%        originalFreq/boundsRefined for the detector's TRUE original
%        bounds and whether refinement actually happened.
%
% B. Miller, AAD, 2026

fprintf('\n=== QC: %s ===\n', cfg.name);

f = fullfile(cfg.folder, cfg.filePattern);
if ~exist(f,'file')
    warning('Detection file not found for %s: %s', cfg.name, f);
    d = table();
    return;
end

wavFolder = '';
if isfield(cfg.siteMeta,'wavFolder')
    wavFolder = cfg.siteMeta.wavFolder;
end
if isempty(wavFolder)
    warning(['No wavFolder in cfg.siteMeta for %s -- soundFolder will be ' ...
        'empty, SNR estimation will fail.'], cfg.name);
end

d = cfg.converter(f, wavFolder, cfg.siteCode);
d.site = repmat(string(cfg.siteCode), height(d), 1);
d.datetime = datetime(d.t0, 'ConvertFrom', 'datenum');
fprintf('%s: %d raw detections\n', cfg.name, height(d));

%% Validate + attach SNR
d = qcValidateDetections(d, cfg.siteMeta);

validRows   = d(d.qcValid, :);
invalidRows = d(~d.qcValid, :);

calib = [];
if useCalibration
    calib = cfg.siteMeta;
end

if ~isempty(validRows)
    % Preserve validRows' OWN bounds (the detector's true original,
    % before any fixedSnrBand/snrTimeBuffer_s widening below) as the
    % originalT0/originalTEnd/originalFreq that ultimately land in d --
    % NOT whatever qcAttachSnr itself sees as "original", since
    % snrInputRows' bounds get overwritten by the widening BEFORE
    % qcAttachSnr ever runs. qcAttachSnr's own originalXxx columns would
    % otherwise just capture the widened search window, not the real
    % detector output.
    trueOriginalT0   = validRows.t0;
    trueOriginalTEnd = validRows.tEnd;
    trueOriginalFreq = validRows.freq;

    snrInputRows = validRows;
    if isfield(cfg,'fixedSnrBand') && ~isempty(cfg.fixedSnrBand)
        % Use a known-physical band for SNR estimation rather than the
        % detector's own (possibly unreliable) per-detection freq bounds
        % -- valid when the call type's frequency content is known a
        % priori. Only affects the SNR call; validRows/d keep the
        % detector-reported fLow/fHigh untouched for QC and diagnostics
        % (restored below from trueOriginalFreq regardless of what
        % qcAttachSnr does to snrInputRows' own copy).
        snrInputRows.fLow  = repmat(cfg.fixedSnrBand(1), height(snrInputRows), 1);
        snrInputRows.fHigh = repmat(cfg.fixedSnrBand(2), height(snrInputRows), 1);
        snrInputRows.freq  = repmat(cfg.fixedSnrBand, height(snrInputRows), 1);
    end
    if isfield(cfg,'snrTimeBuffer_s') && ~isempty(cfg.snrTimeBuffer_s) && cfg.snrTimeBuffer_s > 0
        % Widen the audio window fed into SNR estimation only -- bounds
        % won't be tight (not meant to be), but short/fixed-frame
        % durations otherwise don't leave enough signal for a stable
        % estimate. validRows/d keep the original reported t0/tEnd
        % (restored below).
        bufDays = cfg.snrTimeBuffer_s / 86400;
        snrInputRows.t0   = snrInputRows.t0   - bufDays;
        snrInputRows.tEnd = snrInputRows.tEnd + bufDays;
    end

    snrArgs = {'snrType', cfg.snrType, 'nfft', cfg.nfft};
    if isfield(cfg,'nOverlap') && ~isempty(cfg.nOverlap)
        snrArgs = [snrArgs, {'nOverlap', cfg.nOverlap}];
    end
    if isfield(cfg,'ridgeParams') && ~isempty(cfg.ridgeParams)
        % Passed through as-is to snrEstimate -- e.g. kalmanRidge's
        % globalReferencePercentile/meanAbsentDuration_frames/etc.
        snrArgs = [snrArgs, {'ridgeParams', cfg.ridgeParams}];
    end
    if ~isempty(calib)
        snrArgs = [snrArgs, {'calibration', calib}];
    end
    snrResult = qcAttachSnr(snrInputRows, snrArgs{:});

    validRows.snr = snrResult.snr;
    validRows.RL  = snrResult.RL;
    validRows.NL  = snrResult.NL;

    % Refined bounds (2026-07-24): qcAttachSnr already overwrote
    % snrResult's own t0/tEnd/freq in place (refined where
    % boundsRefined==true, else equal to whatever snrInputRows had --
    % the WIDENED search window, not the true original). Bring the
    % refined values into validRows, but restore the TRUE original
    % bounds captured above, not qcAttachSnr's notion of "original".
    if all(ismember({'t0','tEnd','freq','boundsRefined'}, snrResult.Properties.VariableNames))
        validRows.t0            = snrResult.t0;
        validRows.tEnd          = snrResult.tEnd;
        validRows.freq          = snrResult.freq;
        validRows.boundsRefined = snrResult.boundsRefined;
        validRows.originalT0    = trueOriginalT0;
        validRows.originalTEnd  = trueOriginalTEnd;
        validRows.originalFreq  = trueOriginalFreq;
        % Keep fLow/fHigh in sync with the new freq -- qcSummaryReport
        % and others read these directly, not derived from freq.
        if all(ismember({'fLow','fHigh'}, validRows.Properties.VariableNames))
            validRows.originalFLow = validRows.fLow;   % true original, before overwrite
            validRows.originalFHigh = validRows.fHigh;
            validRows.fLow  = snrResult.freq(:,1);
            validRows.fHigh = snrResult.freq(:,2);
        end
    end
else
    validRows.snr = nan(0,1); validRows.RL = nan(0,1); validRows.NL = nan(0,1);
end

invalidRows.snr = nan(height(invalidRows),1);
invalidRows.RL  = nan(height(invalidRows),1);
invalidRows.NL  = nan(height(invalidRows),1);
% invalidRows never went through SNR estimation -- give them matching
% originalXxx/boundsRefined columns (pass-through, unrefined) so the
% concatenation below has consistent columns across both subsets.
if any(strcmp('boundsRefined', validRows.Properties.VariableNames))
    invalidRows.boundsRefined = false(height(invalidRows),1);
    invalidRows.originalT0    = invalidRows.t0;
    invalidRows.originalTEnd  = invalidRows.tEnd;
    invalidRows.originalFreq  = invalidRows.freq;
    if any(strcmp('originalFLow', validRows.Properties.VariableNames))
        invalidRows.originalFLow  = invalidRows.fLow;
        invalidRows.originalFHigh = invalidRows.fHigh;
    end
end

d = sortrows([validRows; invalidRows], 't0');

%% Prominent sanity check -- SNR yield, min/median/max for all time-freq/SNR fields
qcSummaryReport(d, 'name', cfg.name);

%% Pattern plausibility: detections/day, flagged categories overlaid
qcDetectionRatePlot(d, 'siteCol','site', 'dateCol','datetime', 'flagCol','qcFlag');
sgtitle(sprintf('%s -- detection rate pattern check', cfg.name), 'Interpreter','none');

%% Valid vs flagged distributions
qcFlagDistributionPlots(d);
sgtitle(sprintf('%s -- valid vs flagged distributions', cfg.name), 'Interpreter','none');

%% SNR / NL / RL diagnostics: monthly time series + histogram
qcSnrDiagnosticPlots(d);
sgtitle(sprintf('%s -- SNR diagnostics', cfg.name), 'Interpreter','none');

%% Visual review: full detection set
fprintf('\n--- %s: %d detections to review: Press Esc to continue ---\n', ...
    cfg.name, height(d));
figure('position',[100, 100, 1280, 720]);
qcPageDetections(d, wavFolder, cfg.spectroParams);

end
