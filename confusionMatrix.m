function [cm, truePosRate, falseAlarmRate, duplicateMatches, precision, recall, numGroundTrue] = confusionMatrix(cap, numAbsent)
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

groundTrue = cap.detect_t1; % Logical index of ground truth
detectTrue = cap.detect_t2;

% All unique detections from t1 are considered the ground true presence 
numGroundTrue = length( unique( cap.key(groundTrue) ) ); 
duplicateMatches = length(cap.key) - length(unique(cap.key));

truePositive = sum(groundTrue & detectTrue);     
falsePositive = sum(detectTrue & ~groundTrue);   
falseNegative = numGroundTrue - truePositive - duplicateMatches;
trueNegative = numAbsent - falsePositive;

cm(1,1) = truePositive;     % True positives
cm(1,2) = falsePositive;    % False positives
cm(2,1) = falseNegative;    % Missed detection (false negative)
cm(2,2) = trueNegative;   % True negative

cm = array2table(cm,'VariableNames',{'presence','absence'},...
    'RowNames',{'predicted presence', 'predicted absence'});

truePosRate = (truePositive-duplicateMatches)./numGroundTrue;
falseAlarmRate = falsePositive./numAbsent;
recall = truePosRate;
precision = truePositive./(truePositive + falsePositive);
