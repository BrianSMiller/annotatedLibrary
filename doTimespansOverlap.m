function overlap = doTimespansOverlap(startTime1, endTime1, startTime2, endTime2) 
% Check and see if two time intervals overlap. 
% if (endTime1 < startTime2 | startTime1 > endTime2)
%     overlap = false;
% else
%     overlap = true;
% end
overlap = ~(endTime1 < startTime2 | startTime1 > endTime2);

