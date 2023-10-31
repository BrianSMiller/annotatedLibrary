function roc = rocDetector(cap, threshold, numTests)
% Convert a cell array of capture history tables to an ROC curve (true
% positive rate, false alarm rate, and duplicate matches).
%
% CAP is a cell array of capture history tables (see captureHistoryTable.m)
% Each capture history table in CAP should correspond to a comparison
% between a 'ground truth' and a test table (e.g. manual annotations and an
% automated acoustic detector).

for i = 1:length(cap)
        c = cap{i};
       
        groundTrue = c.detect_t1;
        detectTrue = c.detect_t2;
        
        numGroundTrue(i) = sum(groundTrue);% height(ravenTable)-length(ignoreIx);
        
        tTrue = (duration*numGroundTrue(i));
        tFalse =  tTotal - tTrue;
        numTests(i) = tFalse / duration;
        
        duplicateMatches(i,j) = length(c.key) - length(unique(c.key));
        truePosRate(i,j) = (sum(groundTrue & detectTrue)-duplicateMatches(i,j))./numGroundTrue(i,j);
        falseAlarmRate(i,j) = sum(detectTrue & ~groundTrue)./numTests(i,j);
end
    
[threshold, sortedIx] = sortrows(threshold(:));
truePosRate = truePosRate(sortedIx,:);
falseAlarmRate = falseAlarmRate(sortedIx,:);
duplicateMatches = duplicateMatches(sortedIx,:);


roc = table(falseAlarmRate,truePosRate,threshold, duplicateMatches);