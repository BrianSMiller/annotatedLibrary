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
%   detTable    - input table with snr, RL, NL columns appended. RL/NL
%                 are calibrated (dB re 1 uPa) if calibration was supplied,
%                 otherwise uncalibrated RMS dB.
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

result = snrEstimate(detTable, params);

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
end
