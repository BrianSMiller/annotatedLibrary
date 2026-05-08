function [overlap, pct1, pct2] = timespanOverlap(start1, end1, start2, end2) 
% overlap is the amount of time that overlasp between two time intervals.
% pct1 is the percentage of overlap compared to the first sound;
% pct2 is the perecntage of overlap compared to the second sound;
if ~doTimespansOverlap(start1, end1, start2, end2)
    overlap = zeros(size(start2));
    pct1 = overlap;
    pct2 = overlap;
%     start = start1;
%     stop = end1;
    return;
end

% They overlap. 
% Sort them so that start1 is always first
if (start2 < start1)
    [start2 end2, start1 end1] = swapTimes(start1, end1, start2, end2);
end

% Check if   timespan2 is entirely contained within timespan1
if end1 > end2 
    overlap = (end2 - start2);
else % Partial overlap between the two
    overlap = (end1 - start2);
end

% Only calculate if requested
if nargout > 1
    pct1 = overlap/(end1-start1)*100;
    pct2 = overlap/(end2-start2)*100;
end

function [s2, e2, s1, e1] = swapTimes(s1, e1, s2, e2)
