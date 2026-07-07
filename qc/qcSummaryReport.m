function qcSummaryReport(d, varargin)
% qcSummaryReport  Print a prominent sanity-check summary: how many valid
% detections actually got a usable SNR, plus min/median/max for duration,
% frequency bounds, bandwidth, and SNR/NL/RL. Meant to surface problems
% like "only 14% of valid detections got an SNR value" immediately,
% rather than needing a follow-up manual check that might not happen.
%
%   qcSummaryReport(d)
%   qcSummaryReport(d, 'name','Coni_LFDCS_Dcall', 'snrWarnThreshold', 0.9)
%
% Name-value args
%   name              - label for the header (default '')
%   snrWarnThreshold  - warn (loudly, via warning()) if the fraction of
%                       valid detections with a non-NaN snr falls below
%                       this (default 0.9)
%
% B. Miller, AAD, 2026

p = inputParser;
addParameter(p,'name','');
addParameter(p,'snrWarnThreshold',0.9);
parse(p,varargin{:});
opt = p.Results;

nTotal = height(d);
hasValidCol = any(strcmp('qcValid',d.Properties.VariableNames));
if hasValidCol
    nValid = sum(d.qcValid);
else
    nValid = nTotal;
end

bar = repmat('=',1,60);
fprintf('\n%s\n', bar);
if ~isempty(opt.name)
    fprintf('QC SUMMARY: %s\n', opt.name);
else
    fprintf('QC SUMMARY\n');
end
fprintf('%s\n', bar);
fprintf('Total detections:    %d\n', nTotal);
if hasValidCol
    fprintf('Valid detections:    %d (%.1f%%)\n', nValid, 100*nValid/max(nTotal,1));
end

hasSnr = any(strcmp('snr',d.Properties.VariableNames));
if hasSnr
    snrOk = sum(~isnan(d.snr));
    pctSnrOk = 100 * snrOk / max(nValid,1);
    fprintf('Valid w/ usable SNR: %d of %d (%.1f%% of valid)\n', snrOk, nValid, pctSnrOk);
    if pctSnrOk/100 < opt.snrWarnThreshold
        warning('qcSummaryReport:lowSnrYield', ...
            ['Only %.1f%% of valid detections got a usable SNR estimate ' ...
             '(threshold %.0f%%). Check nfft/nOverlap/snrTimeBuffer_s ' ...
             'against typical detection duration.'], pctSnrOk, opt.snrWarnThreshold*100);
    end
end

fprintf('\n%-10s %10s %10s %10s %10s\n', 'Field','Min','Median','Max','N (finite)');
fprintf('%s\n', repmat('-',1,60));

fields = {'duration','fLow','fHigh','snr','RL','NL'};
for i = 1:numel(fields)
    fn = fields{i};
    if ~any(strcmp(fn, d.Properties.VariableNames))
        continue;
    end
    x = d.(fn);
    xf = x(isfinite(x));
    if isempty(xf)
        fprintf('%-10s %10s %10s %10s %10d\n', fn, 'NaN','NaN','NaN', 0);
    else
        fprintf('%-10s %10.3g %10.3g %10.3g %10d\n', fn, min(xf), median(xf), max(xf), numel(xf));
    end
end

if all(ismember({'fLow','fHigh'}, d.Properties.VariableNames))
    bw = d.fHigh - d.fLow;
    bwf = bw(isfinite(bw));
    if ~isempty(bwf)
        fprintf('%-10s %10.3g %10.3g %10.3g %10d\n', 'bandwidth', min(bwf), median(bwf), max(bwf), numel(bwf));
    end
end

if hasValidCol && any(strcmp('qcFlag',d.Properties.VariableNames)) && any(~d.qcValid)
    fprintf('\nFlag breakdown:\n');
    disp(groupsummary(d(~d.qcValid,:), 'qcFlag'));
end

fprintf('%s\n\n', bar);
end
