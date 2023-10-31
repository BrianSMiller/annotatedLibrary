function s = annotationSummary(a)
%s = annotationSummary(d) - Compute mean, std of parameters of annotations
% d is a Nx1 array of detection data structures (e.g. from
% ravenTableToDetection)
% s is a Matlab structcontaining mean, std, and 5th/95th percentiles of 
% the duration, frequencies, bandwidth of the annotation


s.siteCode = count({a.siteCode});
s.classification = count({a.classification});
s.duration = stats([a.duration]);
s.fLow = stats([a.fLow]);
s.fHigh = stats([a.fHigh]);
bandwidth = [a.fHigh] - [a.fLow];
s.bandwidth = stats(bandwidth);
s.numberOfAnnotations = length(a);

function s = stats(vec)
s.mean = mean(vec);
s.std = std(vec);
p = prctile(vec(:),[5 95]);
s.p05 = p(:,1);
s.p95 = p(:,2);

function c = count(cellVec)
[uniqueVal, ~, ix]=unique(cellVec);
count = histc(ix, 1:numel(uniqueVal));
c = {uniqueVal{:}, count(:)};