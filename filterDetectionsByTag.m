function t = filterDetectionsByTag(t, tagValue)
% filterDetectionsByTag  Restrict a detection table to rows matching a
% single Tags value.
%
% Meant to be composed with an existing converter via an anonymous
% function closure, e.g.:
%   cfg.converter = @(f,wf,sc) filterDetectionsByTag(kooguTableToDetection(f,wf,sc), 'Bm-D');
%
% This lets one raw Koogu output file (with multiple call types mixed
% together via the Tags column) drive multiple qcDetector.m configs, one
% per call type, without modifying qcDetector.m's fixed 3-argument
% converter call signature, and without touching the existing
% kooguTableToDetection.m converter.
%
% INPUTS
%   t         detection table with a Tags column (as produced by
%             kooguTableToDetection or any other converter)
%   tagValue  the Tags value to keep, e.g. 'Bm-D' or 'Bp-Downsweep'
%
% OUTPUT
%   t  same table, restricted to rows where Tags == tagValue
%
% B. Miller, AAD, 2026

if ~any(strcmp('Tags', t.Properties.VariableNames))
    error('filterDetectionsByTag:noTagsColumn', ...
        'No ''Tags'' column found -- nothing to filter on.');
end

mask = strcmp(string(t.Tags), string(tagValue));
fprintf('filterDetectionsByTag: keeping %d of %d rows where Tags==''%s''\n', ...
    sum(mask), height(t), tagValue);
t = t(mask, :);

end
