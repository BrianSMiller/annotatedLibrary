function c = judgeDetections(c,sampleRate,nfft,noverlap,pre,post,startIndex)
% c = judgeAllDetections(c,sampleRate,nfft,noverlap,pre,post)
% Manually review and judge detections from captureHistoryTable, c.
% Audio will be resampled to sampleRate, and spectrograms will be created
% using nfft point FFTs, with noverlap points of overlap between slices.
% Pre and post are the number of seconds of audio to load before and after
% c.t0 and c.tEnd.
cmPerPixel = 0.10;
nDet = height(c);
wavInfo = wavFolderInfo(c.soundFolder_table1{1});


figure('position',[680    500   679   400]);
set(gcf,'units','centimeters');
pos = get(gcf,'pos');

for i = startIndex:nDet
    try
        if c.verdict(i); continue; end
        % calculate audio times
        startTime = min(c.t0_table1(i), c.t0_table2(i)) - pre/86400;
        endTime = max(c.tEnd_table1(i),  c.tEnd_table2(i)) + post/86400;

        wav =  downsampleSingleChannelFromFiles(wavInfo,startTime,endTime,1,sampleRate);
        [s,f,t,p] = spectrogram(wav,nfft,noverlap,[],sampleRate,'yAxis');

        pdB = 20*log10(p);
        fIx = double(f > 30 & f < 125);
        nanIx = find(fIx==0);
        fIx(nanIx)=nan(size(nanIx));
        pIx = repmat(fIx,1,length(t));
        h = pcolor(t,f,pdB.*pIx);
        set(h,'lineStyle','none');
        %         view([0,90]);
        clim = caxis;
        caxis([-50 0]+clim(2));
        %         caxis([100,200]-200);
        %         ylim([0,0.5]);
        ti = sprintf('Row % 4g/% 4g: %s', i, nDet, datestr(c.t0(i)));
        title(ti);
        pos(3) = size(s,2)*cmPerPixel;
        set(gcf,'position',pos);
        drawnow;
        [~,~,button] = ginput(1);
        switch button
            case 1 % Left click for true positive;
                c.verdict(i) = 1;
            otherwise
                c.verdict(i) = 0;
        end
        disp(ti);
    catch
        fprintf('Exiting on index %g\n',i);
        return
    end
end
end

%%
function c = addSnr(c,params)
% Simple wrapper to around annotationSNR that adds required columns to a
% captureHistoryTable
if nargin < 2
    params.freq = [30 115];
    params.noiseDelay = 0;
end
c.duration = max(c.duration_table1,c.duration_table2);
c.tEnd = c.t0 + c.duration/86400;
c.freq = repmat([20,115],height(c),1);
c.soundFolder = repmat(c.soundFolder_table1(1),height(c),1);
[c.snrLurton, c.rmsSignal, c.rmsNoise, c.noiseVar] = annotationSNR(c,params);
c.snrSimple = 20*log10(c.rmsSignal/c.rmsNoise);
end
