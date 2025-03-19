function pageDetections(c,wavFolder,sp)
% c = readtable('captureHistory_casey2019MGA_vs_denseNetBmD24_judgedBSM_cut.csv');
% wavFolder = 'm:\annotatedLibrary_AAD\Casey2019\wav\';
% sp: Spectrogram parameters 
nDet = height(c);
% % wavInfo = wavFolderInfo(c.soundFolder_table1{1});
wavInfo = wavFolderInfo(wavFolder)
% saveFile = 'captureHistory_casey2019MGA_vs_denseNetBmD17_judgedBSM.csv'


%% View spectrogram clips of detections
if nargin < 3
    sp = spectroParams('blue');
end
sampleRateOrig = wavInfo(1).sampleRate;

button = 1;
k = 14000;

while button ~=27 
% for k = 1:nrow*ncol:height(c) % outer loop for pages
    hTile = tiledlayout(sp.nrow,sp.ncol); % inner loop for plots
    
    for i = 1:min(sp.nrow*sp.ncol,nDet-k)
        nexttile
        ix = k+i-1;

        % uncomment to use start and end times
%         startTime = min(c.t0_table1(ix), c.t0_table2(ix)) - pre/86400;
%         endTime = max(c.tEnd_table1(ix),  c.tEnd_table2(ix)) + post/86400;
            
        % use a 10 s window around the mid-point of the call
%         startTime = min(c.t0_table1(ix), c.t0_table2(ix));
%         endTime = max(c.tEnd_table1(ix),  c.tEnd_table2(ix));
        startTime = c.t0(ix);
        endTime = c.tEnd(ix);
        midpoint = mean([startTime endTime]);
        c.FileOffset_s_(ix)  = c.FileOffset_s_(ix) - sp.pre;
        c.DeltaTime_s_(ix) =  c.DeltaTime_s_(ix) + sp.post;

%         wav =  downsampleSingleChannelFromFiles(wavInfo,startTime,endTime,1,sampleRate);
        wav = downsampleFromWav(c(ix,:),sampleRateOrig,sp.sampleRate);
        [s,f,t,p] = spectrogram(wav,sp.nfft,sp.noverlap,sp.nfft,sp.sampleRate,'yAxis');
        
        pdB = 20*log10(p);
        fIx = double(f > sp.lowFreq & f < sp.highFreq);
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
%         title(sprintf('Row % 4g: %s, Manual/Auto %g/%g',...
%             ix, datestr(c.t0(ix)), c.detect_t1(ix),c.detect_t2(ix)));
        title(sprintf('Detection #%g',ix));
%         text(5,2,sprintf('Analyst: %s; DNN: %s; Judge: %s',...
%              label{c.detect_t1(ix)+1},...
%              label{c.detect_t2(ix)+1},...
%              label{c.verdict(ix)+1} ) ,...
%              'horizontalAlign','center','verticalAlign','baseline');
        xlabel(sprintf(datestr(c.t0(ix))));

    end
    [~,~,button] = ginput(1);
    set(hTile.Title,'String','Seeking');
    drawnow;
    switch button
        case 27 % Esc key (quit)
            set(hTile.Title,'String','Done');
        case 28 % Left arrow (previous page)
            k = k-sp.nrow*sp.ncol;
        case 29 % Right arrow (next page)
            k = k+sp.nrow*sp.ncol;
        case 30 % Up arrow (foward a big increment like 5%)
            k = k + floor(nDet*.05);
        case 31 % Down arrow (back a big increment like %5)
            k = k - floor(nDet*.05);    
        case 32 % Spacebar (nothing)
        otherwise % (next page)
            k = k+sp.nrow*sp.ncol;
    end
    % Wrap around if jumped past end or beginning
    k = rem(k,nDet);
    if k < 1
        k = nDet + k;
    end
end

function wav = downsampleFromWav(c,inputRate,outputRate)
startSample = c.FileOffset_s_*inputRate;
endSample = startSample + (c.DeltaTime_s_*inputRate);
wav = audioread(fullfile(c.soundFolder{1},c.BeginFile{1}),floor([startSample, endSample]));
q =  inputRate/outputRate;
wav = resample(wav,1,q);
