function overlap = doTimespansOverlap(startTime1, endTime1, startTime2, endTime2, buffer)
% Check if two time intervals overlap, with an optional time buffer.
% A buffer expands each interval by +/- buffer before checking overlap.
%
% Usage:
%   overlap = doTimespansOverlap(startTime1, endTime1, startTime2, endTime2)
%   overlap = doTimespansOverlap(startTime1, endTime1, startTime2, endTime2, buffer)

if nargin < 5
    buffer = 0;
end

overlap = ~(endTime1 + buffer < startTime2 - buffer | startTime1 - buffer > endTime2 + buffer);

