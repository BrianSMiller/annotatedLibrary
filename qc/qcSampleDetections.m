function sampleTab = qcSampleDetections(detTable, n, varargin)
% qcSampleDetections  Draw a QC sample from a detection table, stratified
% by a chosen variable (default: snr) so that low, mid, and high values
% are all represented -- avoids a random sample skewing toward the
% common (usually mid-range) case, which is exactly the failure mode of
% just eyeballing the top of a sorted table.
%
%   sampleTab = qcSampleDetections(detTable, n)
%   sampleTab = qcSampleDetections(detTable, n, 'stratifyBy','snr', ...
%                   'nBins',5,'sortBy','snr','shuffle',false)
%
% Inputs
%   detTable   - table of detections, one row per detection
%   n          - total number of detections to sample (roughly divided
%                across bins)
%
% Name-value args
%   stratifyBy - column to stratify by (default 'snr'). Falls back to
%                unstratified random sample with a warning if missing.
%   nBins      - number of strata (default 5)
%   sortBy     - column to sort output by (default = stratifyBy).
%                Pass '' to leave sample order as drawn.
%   shuffle    - if true, shuffle rows instead of sorting (default false)
%
% Output
%   sampleTab  - subset of detTable, ordered per sortBy/shuffle. Feed
%                directly into qcPageDetections.
%
% B. Miller, AAD, 2026

p = inputParser;
addParameter(p,'stratifyBy','snr');
addParameter(p,'nBins',5);
addParameter(p,'sortBy',[]);
addParameter(p,'shuffle',false);
parse(p,varargin{:});
opt = p.Results;
if isempty(opt.sortBy)
    opt.sortBy = opt.stratifyBy;
end

nDet = height(detTable);
n = min(n,nDet);

if ~any(strcmp(opt.stratifyBy,detTable.Properties.VariableNames))
    warning('qcSampleDetections:noStratifyColumn', ...
        '"%s" not found in detTable -- falling back to unstratified random sample.', ...
        opt.stratifyBy);
    ix = randperm(nDet,n);
    sampleTab = detTable(ix,:);
else
    edges = quantile(detTable.(opt.stratifyBy), linspace(0,1,opt.nBins+1));
    edges(1) = -inf; edges(end) = inf; % catch full range incl. outliers
    binIx = discretize(detTable.(opt.stratifyBy), edges);
    perBin = ceil(n/opt.nBins);
    keepIx = [];
    for b = 1:opt.nBins
        binRows = find(binIx==b);
        if isempty(binRows); continue; end
        nTake = min(perBin,length(binRows));
        keepIx = [keepIx; binRows(randperm(length(binRows),nTake))]; %#ok<AGROW>
    end
    if length(keepIx) > n
        keepIx = keepIx(randperm(length(keepIx),n));
    end
    sampleTab = detTable(keepIx,:);
end

if opt.shuffle
    sampleTab = sampleTab(randperm(height(sampleTab)),:);
elseif ~isempty(opt.sortBy) && any(strcmp(opt.sortBy,sampleTab.Properties.VariableNames))
    sampleTab = sortrows(sampleTab,opt.sortBy);
end
end
