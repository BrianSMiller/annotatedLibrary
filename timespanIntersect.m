function [start stop] = timespanIntersect(start1, end1, start2, end2) 
% Find the intersection of two timespans.
% start1 and end1 define the start and stop times of the first timespan,
% and are matlab datenums.
% Returns: start and stop of the intersection of two timespans 
start = [];
stop = [];
if ~doTimespansOverlap(start1, end1, start2, end2)
    return;
end

% They overlap. 
% Sort them so that start1 is always first
if (start2 < start1)
    [start2 end2, start1 end1] = swapTimes(start1, end1, start2, end2);
end
start = start2;

% Check if timespan2 is entirely contained within timespan1
if end1 > end2 
    stop = end2;
else % Partial overlap between the two
    stop = end1;
end
if stop < start
    disp('warning:')
    disp(datestr([start;stop]));
    [stop, start] = swapStartStop(start,stop);
end

function [s2, e2, s1, e1] = swapTimes(s1, e1, s2, e2)

function [e, s] = swapStartStop(s, e)
