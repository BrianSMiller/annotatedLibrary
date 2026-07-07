function qcPrintDetections(detTable, varargin)
% qcPrintDetections  Print detections with human-readable datetimes for
% quick console inspection -- e.g. a handful of rows from one flag
% category, to look at actual values rather than aggregate counts or
% histograms when a diagnostic plot doesn't tell the full story.
%
%   qcPrintDetections(detTable)
%   qcPrintDetections(detTable, 'n', 10, 'flagValue', 'OutsideDeploymentWindow')
%   qcPrintDetections(d(d.qcFlag=="NonIncreasingFreq",:), 'n', 20)
%
% Inputs
%   detTable   - detection table (any subset -- filter yourself, or use
%                flagValue below)
%
% Name-value args
%   n          - max rows to print (default 10)
%   flagValue  - if set, filters to detTable.(flagCol) == flagValue first
%   flagCol    - column to filter on (default 'qcFlag')
%
% B. Miller, AAD, 2026

p = inputParser;
addParameter(p,'n',10);
addParameter(p,'flagValue','');
addParameter(p,'flagCol','qcFlag');
parse(p,varargin{:});
opt = p.Results;

if ~isempty(opt.flagValue)
    if ~any(strcmp(opt.flagCol,detTable.Properties.VariableNames))
        error('qcPrintDetections:noFlagCol', '"%s" not found in detTable.', opt.flagCol);
    end
    detTable = detTable(string(detTable.(opt.flagCol)) == string(opt.flagValue), :);
end

if isempty(detTable)
    fprintf('qcPrintDetections: no rows match.\n');
    return;
end

n = min(opt.n, height(detTable));
sub = detTable(1:n, :);

printCols = intersect({'site','siteCode','t0','tEnd','duration','fLow','fHigh', ...
    'snr','RL','NL','qcFlag'}, sub.Properties.VariableNames, 'stable');

out = table();
for c = 1:numel(printCols)
    col = printCols{c};
    if any(strcmp(col,{'t0','tEnd'}))
        out.(col) = cellstr(datestr(sub.(col), 'yyyy-mm-dd HH:MM:SS.FFF'));
    else
        out.(col) = sub.(col);
    end
end

fprintf('qcPrintDetections: showing %d of %d row(s)\n', n, height(detTable));
disp(out);
end
