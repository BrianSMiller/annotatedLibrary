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
%        with qcValid/qcFlag columns from qcValidateDetections
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
    snrInputRows = validRows;
    if isfield(cfg,'fixedSnrBand') && ~isempty(cfg.fixedSnrBand)
        % Use a known-physical band for SNR estimation rather than the
        % detector's own (possibly unreliable) per-detection freq bounds
        % -- valid when the call type's frequency content is known a
        % priori. Only affects the SNR call; validRows/d keep the
        % detector-reported fLow/fHigh untouched for QC and diagnostics.
        snrInputRows.fLow  = repmat(cfg.fixedSnrBand(1), height(snrInputRows), 1);
        snrInputRows.fHigh = repmat(cfg.fixedSnrBand(2), height(snrInputRows), 1);
        snrInputRows.freq  = repmat(cfg.fixedSnrBand, height(snrInputRows), 1);
    end
    if isfield(cfg,'snrTimeBuffer_s') && ~isempty(cfg.snrTimeBuffer_s) && cfg.snrTimeBuffer_s > 0
        % Widen the audio window fed into SNR estimation only -- bounds
        % won't be tight (not meant to be), but short/fixed-frame
        % durations otherwise don't leave enough signal for a stable
        % estimate. validRows/d keep the original reported t0/tEnd.
        bufDays = cfg.snrTimeBuffer_s / 86400;
        snrInputRows.t0   = snrInputRows.t0   - bufDays;
        snrInputRows.tEnd = snrInputRows.tEnd + bufDays;
    end

    snrArgs = {'snrType', cfg.snrType, 'nfft', cfg.nfft};
    if isfield(cfg,'nOverlap') && ~isempty(cfg.nOverlap)
        snrArgs = [snrArgs, {'nOverlap', cfg.nOverlap}];
    end
    if ~isempty(calib)
        snrArgs = [snrArgs, {'calibration', calib}];
    end
    snrResult = qcAttachSnr(snrInputRows, snrArgs{:});

    validRows.snr = snrResult.snr;
    validRows.RL  = snrResult.RL;
    validRows.NL  = snrResult.NL;
else
    validRows.snr = nan(0,1); validRows.RL = nan(0,1); validRows.NL = nan(0,1);
end

invalidRows.snr = nan(height(invalidRows),1);
invalidRows.RL  = nan(height(invalidRows),1);
invalidRows.NL  = nan(height(invalidRows),1);

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
