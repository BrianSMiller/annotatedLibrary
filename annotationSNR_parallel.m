% function result =  annotationSNR_parallel(annot,params)
% % Run annotationSNR in parallel using the parallel processing toolbox. 
% % TODO: Check for required parallel processing toolbox and other 
% % dependencies
% 
% if istable(annot)
%     annot = table2struct(annot);
% end
% nDet = length(annot);
% 
% if isempty(gcp('nocreate'))
%     parpool('Processes', feature('numcores')-1);
% end
% 
% result = nan(nDet,1); rmsSignal=result; rmsNoise=result; noiseVar=result;
% parfor i = 1:nDet
%     [snr(i), rmsSignal(i), rmsNoise(i), noiseVar(i)] = ...
%         annotationSNR(annot(i),params);
%     if rem(i,floor(nDet/100))==0
% 
%     end
% end
% 
% snr = snr(:);
% signalRMSdB = 10*log10(rmsSignal(:));
% noiseRMSdB = 10*log10(rmsNoise(:));
% noiseDev = noiseVar(:);
% result = table(snr,signalRMSdB,noiseRMSdB,noiseDev);
% 
function result = annotationSNR_parallel(annot,params)
% Run annotationSNR in parallel using the parallel processing toolbox.
% TODO: Check for required parallel processing toolbox and other
% dependencies
if istable(annot)
    annot = table2struct(annot);
end
nDet = length(annot);

if isempty(gcp('nocreate'))
    parpool('Processes', feature('numcores')-1);
end

fprintf('SNR Analysis Started at: %s\n', datestr(now))
fprintf('%g annotations to process\n', nDet);
fprintf('Progress (percent)\n');
fprintf('0|-----------|----------50-----------|---------100|\n ');

% Create DataQueue and prepare progress display
D = parallel.pool.DataQueue;
afterEach(D, @(~) fprintf('#'));

result = nan(nDet,1); 
rmsSignal = result; 
rmsNoise = result; 
noiseVar = result;

% Calculate progress increment (every ~2%)
progInc = max(1, floor(nDet/50));

parfor i = 1:nDet
    [snr(i), rmsSignal(i), rmsNoise(i), noiseVar(i)] = ...
        annotationSNR(annot(i), params);
    
    if rem(i, progInc) == 0 || i == nDet
        send(D, i); % Console progress
    end
end

fprintf('\n'); % New line after progress bar

snr = snr(:);
signalRMSdB = 10*log10(rmsSignal(:));
noiseRMSdB = 10*log10(rmsNoise(:));
noiseDev = noiseVar(:);
result = table(snr, signalRMSdB, noiseRMSdB, noiseDev);

fprintf('Completed at: %s\n', datestr(now))