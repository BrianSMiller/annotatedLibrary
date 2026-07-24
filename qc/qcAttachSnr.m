function detTable = qcAttachSnr(detTable, varargin)
% qcAttachSnr  Attach snr/RL/NL columns to a detection table using bsnr's
% snrEstimate(), calibrated if a calibration struct is supplied.
%
% Assumes detTable has already been validated (see qcValidateDetections)
% -- malformed frequency ranges will crash snrEstimate's parfor batch,
% and validity checking is that function's job, not this one's.
%
%   detTable = qcAttachSnr(detTable)
%   detTable = qcAttachSnr(detTable, 'snrType','spectrogramSlices', ...
%                   'calibration', metaDataCasey2019(), 'nfft', 512)
%
% Inputs
%   detTable    - detection table with the fields snrEstimate expects:
%                 soundFolder, t0, tEnd, channel, and freq (Nx2 [fLow fHigh]
%                 per row). If detTable has separate fLow/fHigh columns
%                 instead of a combined freq column, they're combined
%                 automatically.
%
% Name-value args (passed through to snrEstimate; a few pulled out
% explicitly since batch QC wants them fixed across all rows)
%   snrType      - SNR method (default 'spectrogramSlices')
%   calibration  - instrument calibration struct (e.g. metaDataCasey2019());
%                  omit for uncalibrated RMS levels
%   nfft         - fixed FFT length. Always set this explicitly for batch
%                  QC -- bsnr derives nfft per-row from annotation duration
%                  otherwise, which makes SNR non-comparable across rows.
%   ...          - any other snrEstimate params, passed through as-is
%
% Output
%   detTable    - input table with the following appended/updated:
%                   snr, RL, NL  -- as before. RL/NL are calibrated
%                                (dB re 1 uPa) if calibration was
%                                supplied, otherwise uncalibrated RMS dB.
%                   t0, tEnd, freq  -- OVERWRITTEN IN PLACE (2026-07-24)
%                                with bsnr's recovered bounds, for
%                                methods with a bounds-refinement concept
%                                (momentumRidge, kalmanRidge). Unchanged
%                                for every other method, or when
%                                refinement found no call at trackable
%                                confidence (dropout, no fallback) --
%                                deliberately in place, not a separate
%                                column, so every EXISTING downstream
%                                consumer that already reads t0/tEnd/freq
%                                as canonical (qcPageDetections,
%                                qcValidateDetections,
%                                qcFlagDistributionPlots, matchbox, ...)
%                                benefits automatically with no changes
%                                needed there.
%                   originalT0, originalTEnd, originalFreq  -- NEW. The
%                                detector's TRUE original bounds,
%                                preserved before any overwrite above.
%                                Guarded against qcAttachSnr running
%                                twice on the same table -- only
%                                populated if not already present, so a
%                                second call can't clobber the real
%                                original with an already-refined value.
%                   boundsRefined -- NEW. true where t0/tEnd/freq above
%                                are genuinely refined, false where
%                                they're unchanged from originalT0/
%                                originalTEnd/originalFreq.
%
% NOTE: bsnr's documented output field for calibrated signal level is
% signalBandLevel_dBuPa. The doc doesn't explicitly name the calibrated
% noise-level field -- this assumes noiseBandLevel_dBuPa by symmetry.
% Run once on a small table and check fieldnames(result) to confirm; tell
% me the actual name if it differs and I'll fix the one line below.
%
% B. Miller, AAD, 2026

p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'snrType','spectrogramSlices');
addParameter(p,'calibration',[]);
addParameter(p,'nfft',[]);
parse(p,varargin{:});
opt = p.Results;

if ~any(strcmp('freq',detTable.Properties.VariableNames))
    if all(ismember({'fLow','fHigh'},detTable.Properties.VariableNames))
        detTable.freq = [detTable.fLow, detTable.fHigh];
    else
        error('qcAttachSnr:noFreq', ...
            'detTable needs a freq column ([fLow fHigh] per row), or fLow/fHigh columns to build one from.');
    end
end

params = p.Unmatched;
params.snrType = opt.snrType;
if ~isempty(opt.calibration)
    params.calibration = opt.calibration;
end
if ~isempty(opt.nfft)
    params.nfft = opt.nfft;
elseif ~isfield(params,'nfft')
    warning('qcAttachSnr:nfftNotSet', ...
        ['nfft not set -- bsnr will derive it per-row from annotation ' ...
         'duration. Set nfft explicitly for QC batches so SNR stays ' ...
         'comparable across detections.']);
end

result = snrEstimateQuiet(detTable, params);

detTable.snr = result.snr;

if any(strcmp('signalBandLevel_dBuPa',result.Properties.VariableNames))
    detTable.RL = result.signalBandLevel_dBuPa;
else
    detTable.RL = result.signalRMSdB;
end

if any(strcmp('noiseBandLevel_dBuPa',result.Properties.VariableNames))
    detTable.NL = result.noiseBandLevel_dBuPa;
elseif any(strcmp('noiseRMSdB',result.Properties.VariableNames))
    detTable.NL = result.noiseRMSdB;
end

% NEW (2026-07-24): recovered bounds. Preserve the detector's TRUE
% original bounds first (guarded against qcAttachSnr running twice on
% the same table -- only copy if not already present, so a second call
% can't overwrite the real original with an already-refined value), then
% overwrite t0/tEnd/freq in place with bsnr's result. Safe for
% non-refined rows too: result.t0/tEnd/freq already equal the original
% values there (boundsRefined==false), so the overwrite is a no-op.
% Deliberately IN PLACE, not a separate refinedXxx column -- every
% existing downstream consumer (qcPageDetections, qcValidateDetections,
% qcFlagDistributionPlots, matchbox, ...) already reads t0/tEnd/freq as
% canonical; a differently-named column would need each of those
% individually updated to know to prefer it.
if all(ismember({'t0','tEnd','freq','boundsRefined'}, result.Properties.VariableNames))
    if ~any(strcmp('originalT0', detTable.Properties.VariableNames))
        detTable.originalT0   = detTable.t0;
        detTable.originalTEnd = detTable.tEnd;
        detTable.originalFreq = detTable.freq;
    end
    detTable.t0            = result.t0;
    detTable.tEnd           = result.tEnd;
    detTable.freq           = result.freq;
    detTable.boundsRefined  = result.boundsRefined;
    % Keep fLow/fHigh columns (if present) in sync with the new freq --
    % several QC functions (qcSummaryReport in particular) read fLow/
    % fHigh directly rather than deriving them from freq at read-time.
    if all(ismember({'fLow','fHigh'}, detTable.Properties.VariableNames))
        detTable.fLow  = result.freq(:,1);
        detTable.fHigh = result.freq(:,2);
    end
    nRefined = sum(result.boundsRefined);
    fprintf('qcAttachSnr: %d of %d rows got refined bounds (snrType=''%s'')\n', ...
        nRefined, height(detTable), opt.snrType);
end

end

function result = snrEstimateQuiet(detTable, params)
% Suppress kalmanRidge's dropout warnings for the duration of this batch
% call -- both on the client AND on any existing parallel pool workers
% (warning state is per-session; a plain warning('off',...) in the client
% never reaches workers that already existed before this call, which is
% the actual cause if warnings still show up despite this wrapper --
% see this function specifically, not a client-only warning toggle).
warnIds = {'snrKalmanRidge:dropout', 'deriveBoundsFromPresenceTrack:neverPresent'};

pool = gcp('nocreate');
setWorkerWarnings = @(state) cellfun(@(id) warning(state, id), warnIds);

if ~isempty(pool)
    parfevalOnAll(pool, setWorkerWarnings, 0, 'off');
end
setWorkerWarnings('off');   % client session too

result = snrEstimate(detTable, params);

setWorkerWarnings('on');
if ~isempty(pool)
    parfevalOnAll(pool, setWorkerWarnings, 0, 'on');
end

end
