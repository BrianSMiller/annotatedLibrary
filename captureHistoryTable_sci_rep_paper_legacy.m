function [cap, cm] = captureHistoryTable(table1, table2, varargin )
% function [cap, cm] = captureHistoryTable(table1, table2, varargin )
% Create a Capture History Table, CAP, from two tables of detections,
% table1 & table2. Rows in the capture history table represent detections,
% and detection time-frequency bounds from the two tables are compared to
% determine if a detection from table2 overlaps with any from Table1.
% Detections from table2 that match with a detection from Table 1 are then
% combined into a single row within CAP. CaptureHistoryTables are often used in
% double platform and mark-recapture analyses.
% 
% Detections tables must have one row for each detection, and must have
% columns called, t0, tEnd, fLow, and fHigh, which represent the start
% time, end time, lower frequency bound, and upper frequency bound
% respectively.
% 
% The columns in the capture history table correspond to those from the
% original detector. Thus, the captureHistoryTable will contain all of the
% columns from both tables for each row. Columns with duplicate names in
% both tables (e.g. t0, tEnd, fLow, fHigh, and others) will have either _t1
% or _t2 appended to their name depending whether they correspond to
% columns from table1 or table2 respectively.
%
% Additional columns, detect_t1, detect_t2, and key are added to CAP by
% this function. detect_t1 and detect_t2 contain logical values depending
% on whether that row contains a detection from table1 or table2
% respectively. Key is an identifier used to match detections between
% tables. Duplicate values of key indicate duplicate matches between table
% 2 and table 1.


% Start the capture history table by adding all of the detections from
% table 1. 
%
% Since tables may have duplicate columns, append the table
% name to the end of each column
% Commented out because Matlab already does this when joining tables.
% t1names = strcat(fieldnames(table1),'_t1')
% t2names = strcat(fieldnames(table2),'_t2')
tic;
str='';
table1.key = [1:height(table1)]';
table2.key = nan(height(table2),1);
table1.detect_t1 = true(size(table1.t0));
table2.detect_t2 = true(size(table2.t0));
table2.overlap = nan(height(table2),1);

% In Matlab 2014a, accessing table cells by index is incredibly slow if
% done within a for loop, so we extract these columns into separate
% variables to speed things up.
start1 = table2.t0; end1 = table2.tEnd;
start2 = table1.t0; end2 = table1.tEnd;
lowFreq1 = table2.fLow; hiFreq1 = table2.fHigh;
lowFreq2 = table1.fLow; hiFreq2 = table1.fHigh;

% Loop over second table, and determine whether any of the
% detections match those of the first table.
% If they match, then add a matching value to the key of table 2,
% otherwise add a new key to table 2 to identify this as a new 'event' that 
% doesn't match with table 1.
k = height(table1);


for j = 1:height(table2);
    s1 = start1(j); e1 = end1(j);
   
    overlapIx = find(doTimespansOverlap(s1, e1, start2, end2));
    overlap = 0; 
    otherMatches = 0;
    
    % Loop over all matches to calculate the overlap in time & frequency
    for i = overlapIx; 1:height(table1);
        if isempty(i);
            continue;
        end
        s2 = start2(i); e2 = end2(i);

       tOverlap(i) = timespanOverlap(s1, e1, s2, e2); % Time overplap (days)
       
       fOverlap(i) = 0;
       if tOverlap(i) > 0; % Shouldn't this always be > 0 since these are only overlapping calls...?
           tOverlap = tOverlap * 86400; % Convert to s;
           l1 = lowFreq1(j); u1 = hiFreq1(j);
           l2 = lowFreq2(i); u2 = hiFreq2(i);
           fOverlap(i) = timespanOverlap(l1,u1,l2,u2); % Frequency overlap (Hz)
           overlap(i) = tOverlap(i) .* fOverlap(i); 
       end
    end
    
    % Check for more than one matching detection in table 1
    ix = 1;
    if length(overlap) >= 1 
        [overlap, ix] = max(overlap);
        % Check if another detection from table 2 already matches
        otherMatches = find(ix==table2.key);
    end
    table2.overlap(j) = overlap;
    

    if length(otherMatches) > 0        
%         warning('Multiple detections from table 2 matched the same detection in table 1');
    end
    
    if overlap > 0 % Match, so add this to an existing row
        table2.key(j) = ix;
    else % Add a new row to our capture history
        k = k+1;
        table2.key(j) = k;
    end
    if rem(j,500)==0 || j == height(table2)
        fprintf(repmat('\b',1,length(str)));
        str = sprintf('%g/%g rows from table2 compared with table 1 in %g s\n',...
                j, height(table2), toc);  
        fprintf(str);
    end
end

% Merge the two tables into a single table SQL-style using an 'outer-join'
% This table will have a row for every row in table 1, and will put
% matching detections from table 2 into those rows. Rows from table 1 will
% be duplicated if multiple rows from table 2 match them.
cap = outerjoin(table1,table2,'keys',{'key','key'},'MergeKeys',true);

fprintf(['%g events\n',...
         '%g table1\n',...
         '%g table2\n',...
         '%g both\n',...
         '%g duplicate matches from table2 to table1\n'],...
        height(cap),...
        height(table1),...
        height(table2),...
        height(table1)+height(table2)-height(cap),...
        height(cap)-length(unique(cap.key)));