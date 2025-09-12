function [precision, recall, numGroundTrue, duplicateMatches] = det_pr(groundTrue, predictTrue, uid)
% [cm] = confusionMatrix(c)
% Given a captureHistoryTable, cap, from an acoustic detector calculate a
% confusion matrix, cm. We assume that detector1 (_t1 from cap) corresponds
% to the true presence, while detector2 (_t2 from cap) corresponds to the
% predicted presence or absence.
%
% CaptureHistoryTables from passive acoustic detectors usually only contain
% information on actual presence (ground truth), but do not include
% information on actual absence (ground false?). Thus the number of true
% negative detections, numAbsent, must be provided separately from cap.

% All unique detections from t1 are considered the ground true presence 
numGroundTrue = length( unique( uid(groundTrue) ) ); 
duplicateMatches = length(uid) - length(unique(uid));

truePositive = sum(groundTrue & predictTrue);     
falsePositive = sum(predictTrue & ~groundTrue);   
falseNegative = numGroundTrue - truePositive - duplicateMatches;

truePosRate = (truePositive-duplicateMatches)./numGroundTrue;
recall = truePosRate;
precision = truePositive./(truePositive + falsePositive);
