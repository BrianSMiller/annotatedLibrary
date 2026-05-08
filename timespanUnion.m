function [start stop] = timespanUnion(start1, end1, start2, end2) 
% Find the union of two timespans.
% start1, end1, and start2,end2 define the start and stop times of two
% spans and are all matlab datenums.
% Returns: start and stop of the union of two timespans 
% The union of two timepans (only if they overlap) 
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
start = start1;

% Check if timespan2 is entirely contained within timespan1
if end1 > end2 
    stop = end1;
else % Partial overlap between the two
    stop = end2;
end
if stop < start
    disp('warning:')
    disp(datestr([start;stop]));
end

function [s2, e2, s1, e1] = swapTimes(s1, e1, s2, e2)
