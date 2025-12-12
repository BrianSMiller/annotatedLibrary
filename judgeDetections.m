function c = judgeDetections(c,sampleRate,nfft,noverlap,pre,post,startIndex, increment, suppressClicks)
% c = judgeAllDetections(c,sampleRate,nfft,noverlap,pre,post)
% Manually review and judge detections from captureHistoryTable, c.
% Audio will be resampled to sampleRate, and spectrograms will be created
% using nfft point FFTs, with noverlap points of overlap between slices.
% Pre and post are the number of seconds of audio to load before and after
% c.t0 and c.tEnd.

%% Create verdict and judged columns unless they already exist
if ~any(strcmp('verdict',fieldnames(c)))
    c.verdict = nan(size(c.detect_observer1));
end
if ~any(strcmp('judged',fieldnames(c)))
    c.judged = false(size(c.verdict));
end

% Program flow options
showJudged = false;     % Set to false to bypass already judged detections
overwriteJudged = false; % subordinate to showJudged; Set false to view-only 
                        % set true to redo-adjudication 

% Visualisation options
cmPerPixel = 0.05;      % Horizontal pixel scaling (cm per pixel).
showObsBoxes = true;    % Show observer detection boundaries
windowDur = 60;         % Duration (s) of audio clip (if using fixed duration)
pos = [100,400,679,400];% Starting position of figure window (pixels)
freqBand = [17 30];     % Lower and upper limits of frequency band used for color scale
cAdjust = [5 99];       % Color scale floor and ceiling percentiles
windowType = 'fixed';
if nargin < 9
    suppressClicks = [];
end

nDet = height(c);
wavInfo = wavFolderInfo(c.soundFolder_observer1{1});

figure('position',pos);
colormap(flipud(gray));
set(gcf,'units','centimeters');

i = startIndex;
while (i <= nDet)
    pos = get(gcf,'pos');
    msg = '';
    try
        if c.judged(i)==true & ~isnan(c.verdict(i))
            if showJudged
                msg = sprintf('Already judged as %g', c.verdict(i));
            else
                i = i+increment; % Uncomment to skip already judged
                continue
            end
        end % keep existing verdicts

        nObservers = sum(cellfun(@isempty,...
            strfind(c.Properties.VariableNames,'detect_observer') ...
            )==0);


        % calculate audio times
        callStart = min(c{i,append("t0_observer",string(1:nObservers))});
        callEnd = max(c{i,append("tEnd_observer",string(1:nObservers))});

        switch(windowType)
            case 'prePost'
                % Fixed pre and post
                startTime = callStart - pre/86400;
                endTime = callEnd + post/86400;
            case 'fixed'
                startTime = mean([callStart callEnd])-0.5*windowDur/86400;
                endTime = mean([callStart callEnd])+0.5*windowDur/86400;
            otherwise
                startTime = callStart;
                endTime = callEnd;
        end
        wav =  downsampleSingleChannelFromFiles(wavInfo,startTime,endTime,1,sampleRate);

        % Remove clicks
        if ~isempty(suppressClicks)
            threshold = 3;
            power = 1000;
            if isfield(suppressClicks,'threshold')
                threshold = suppressClicks.threshold;
            end
            if isfield(suppressClicks,'power')
                power = suppressClicks.power;
            end
            wav = removeClicks(wav,threshold,power);
        end



        [s,f,t,p] = spectrogram(wav,nfft,noverlap,[],sampleRate,'yAxis');

        pdB = 20*log10(p);
        fIx = double(f > freqBand(1) & f < freqBand(2));
        nanIx = find(fIx==0);
        fIx(nanIx)=nan(size(nanIx));
        pIx = repmat(fIx,1,length(t));
        %         h = pcolor(t,f,pdB.*pIx);
        h = pcolor(t,f,pdB);

        set(h,'lineStyle','none');
        %         view([0,90]);
        %         clim(1) = min(pdB(:).*pIx(:));
        clim(1) = prctile(pdB(:).*pIx(:),cAdjust(1));
        %         clim(2) = max(pdB(:).*pIx(:));
        clim(2) = prctile(pdB(:).*pIx(:),cAdjust(2));

        set(gca,'clim',clim);

        %         caxis([100,200]-200);
        %         ylim([0,0.5]);
        horzline([29,24],'c--');
        vertline([(callStart-startTime)*86400 (callEnd-startTime)*86400],'c--')


        if showObsBoxes
            l = c{i,append("t0_observer",string(1:nObservers))}-startTime;
            r = c{i,append("tEnd_observer",string(1:nObservers))}-startTime;
            top = c{i,append("freq_observer",string(1:nObservers),'_1')};
            bot = c{i,append("freq_observer",string(1:nObservers),'_2')};

            ho = ishold;
            hold on;
            obsColors = lines(nObservers);
            for j = 1:nObservers
                plotBox([l(j),r(j)]*86400,[bot(j),top(j)],obsColors(j,:),'linewidth',1.5);
                text(l(j)*86400,top(j),num2str(j),'color',obsColors(j,:))
            end
            if ho==0; hold off; end
        end

        ti = sprintf('#% 4g/%g floor/ceil percentile=[%g %g]\t(Row % 4g: %s)',...
            1+(i-1)/increment,floor(nDet/increment), cAdjust, ...
            i, datestr(c.t0(i)));

        title(ti);
        text(mean(xlim),max(ylim),msg,'VerticalAlignment','top',...
            'FontSize',14,'FontWeight','bold');

        pos(3) = size(s,2)*cmPerPixel;
        set(gcf,'position',pos);
        drawnow;

        [~,~,button] = ginput(1);
        switch button
            case 28 % Left arrow shift color floor down (darker)
                cAdjust = cAdjust+[-1 0];
            case 29
                cAdjust = cAdjust+[1 0]; % Down arrow color floor down (ligher)
            case 30 % Up arrow shift color ceiling up (lighter)
                cAdjust = cAdjust+[0 1];
            case 31
                cAdjust = cAdjust+[0 -1]; % Down arrow to shift cieling down (darker)
            case 1 % Left click for true positive;
                % Only write if new verdicts (unless overwrite is true)
                if c.judged(i)~=true || overwriteJudged
                    c.verdict(i) = 1;
                    c.judged(i)=true;
                end
                i = i + increment;
            case 3 % Right click for false positive
                % Only write new verdicts or if overwriteJudged==true)
                if c.judged(i)~=true || overwriteJudged
                    c.verdict(i) = 0;
                    c.judged(i)=true;
                end
                i = i + increment;
        end
        % Ensure percentile stay within 0 and 100;
        cAdjust(1) = max(cAdjust(1),0);
        cAdjust(2) = min(cAdjust(2),100);

        disp(ti);
    catch
        l = lasterror;
        disp(l);
        for j = 1:length(l.stack); disp(l.stack(j)); end
        fprintf('Exiting on index %g\n',i);
        return
    end
end
