function qcFlagDistributionPlots(detTable, varargin)
% qcFlagDistributionPlots  Compare duration/fLow/fHigh/bandwidth
% distributions between valid and flagged detections. Answers: do
% flagged detections look like a genuinely distinct population
% (systematic detector/export problem), or just the tail of the same
% distribution (occasional outliers not worth much concern)?
%
%   qcFlagDistributionPlots(detTable)
%   qcFlagDistributionPlots(detTable, 'validCol','qcValid', 'nBins', 80)
%
% Inputs
%   detTable  - detection table with duration (or t0/tEnd), fLow/fHigh
%               (or freq), and a logical valid column (default 'qcValid',
%               as produced by qcValidateDetections)
%
% Name-value args
%   validCol  - logical column marking valid detections (default 'qcValid')
%   nBins     - number of histogram bins (default 80)
%
% Uses the table's own 'duration' column if present, rather than
% recomputing from t0/tEnd -- qcValidateDetections may have clamped
% t0/tEnd for borderline-valid rows, and the diagnostic here is about
% what the detector originally reported, not the post-clamp values.
%
% Bandwidth is computed as fHigh-fLow (not any detector-native bandwidth
% column, e.g. LFDCS's Bandwidth_Hz) so this stays generic across detectors.
%
% B. Miller, AAD, 2026

p = inputParser;
addParameter(p,'validCol','qcValid');
addParameter(p,'nBins',80);
parse(p,varargin{:});
opt = p.Results;

if ~any(strcmp(opt.validCol, detTable.Properties.VariableNames))
    error('qcFlagDistributionPlots:noValidCol', ...
        '"%s" not found -- run qcValidateDetections first.', opt.validCol);
end

if ~any(strcmp('fLow',detTable.Properties.VariableNames)) && ...
        any(strcmp('freq',detTable.Properties.VariableNames))
    detTable.fLow  = detTable.freq(:,1);
    detTable.fHigh = detTable.freq(:,2);
end

if any(strcmp('duration',detTable.Properties.VariableNames))
    durSec = detTable.duration;
else
    durSec = (detTable.tEnd - detTable.t0) * 86400;
end

bandwidth = detTable.fHigh - detTable.fLow;

valid = logical(detTable.(opt.validCol));
nValid = sum(valid);
nFlagged = sum(~valid);

panels = {durSec, detTable.fLow, detTable.fHigh, bandwidth};
titles = {'Duration (s)', 'fLow (Hz)', 'fHigh (Hz)', 'Bandwidth (Hz)'};

figure;
tl = tiledlayout(numel(panels),1);
for i = 1:numel(panels)
    nexttile;
    x = panels{i};
    finiteX = x(isfinite(x));
    if isempty(finiteX)
        title(sprintf('%s (no finite values)', titles{i}));
        continue;
    end
    edges = linspace(min(finiteX), max(finiteX), opt.nBins+1);

    histogram(x(valid), edges, 'FaceColor',[0.2 0.4 0.8], 'FaceAlpha',0.6, ...
        'DisplayName', sprintf('Valid (n=%d)', nValid));
    hold on;
    if nFlagged > 0
        histogram(x(~valid), edges, 'FaceColor','r', 'FaceAlpha',0.6, ...
            'DisplayName', sprintf('Flagged (n=%d)', nFlagged));
    end
    hold off;

    title(titles{i});
    xlabel(titles{i});
    ylabel('Count');
    legend('Location','best');
    grid on;
end
tl.Title.String = 'Valid vs flagged detection distributions';
tl.Title.Interpreter='none';
end
