function qcPageDetections(c, wavFolder, sp, k)
% qcPageDetections  Page through spectrogram clips of detections for
% visual QC -- display only, no verdict recording (that's a separate
% concern; see judgeDetections.m for judged-CSV adjudication workflows).
%
%   qcPageDetections(c, wavFolder)
%   qcPageDetections(c, wavFolder, sp, k)
%
% Inputs
%   c         - detection table with t0, tEnd (datenum). If fLow/fHigh
%               (or freq) are present, a box is drawn around the
%               detection's time-frequency bounds. If snr, RL, and/or NL
%               columns are present they are overlaid as text for the
%               "does this match my ear" check.
%               Sample/sort c beforehand with qcSampleDetections.
%   wavFolder - path to the wav folder for this deployment
%   sp        - spectroParams struct (default: spectroParams('fin')).
%               Two extra optional fields, not in the base spectroParams
%               presets, borrowed from judgeDetections.m:
%                 windowType - 'prePost' (default, uses sp.pre/sp.post
%                              either side of t0/tEnd) or 'fixed' (a
%                              constant-width window centred on the
%                              detection midpoint -- useful for comparing
%                              many differently-sized detections at a
%                              consistent visual scale on one page)
%                 windowDur  - window width in seconds, used only when
%                              windowType is 'fixed'
%   k         - starting row index (default 1)
%
% Navigation: Left/Right arrows = page back/forward, Up/Down = jump 5%,
% Space = redraw current page, Esc = quit.
% Info: I = print every detection on the current page to console
%       (reuses qcPrintDetections); doesn't advance the page.
% S     = save the current page as a PNG via saveAsPng; doesn't advance
%         the page. Filename: qcPageDetections_<startRow>.png
%
% B. Miller, AAD, 2026 -- refactored from pageDetections.m; box overlay,
% windowDur, console print, and screenshot borrowed from judgeDetections.m

if nargin < 3 || isempty(sp)
    sp = spectroParams('fin');
end
if nargin < 4 || isempty(k)
    k = 1;
end
if ~isfield(sp,'windowType') || isempty(sp.windowType)
    sp.windowType = 'prePost';
end

nDet = height(c);
wavInfo = wavFolderInfo(wavFolder);

hasSnr  = any(strcmp('snr',c.Properties.VariableNames));
hasRL   = any(strcmp('RL', c.Properties.VariableNames));
hasNL   = any(strcmp('NL', c.Properties.VariableNames));
hasFreq = all(ismember({'fLow','fHigh'}, c.Properties.VariableNames)) || ...
          any(strcmp('freq', c.Properties.VariableNames));
if any(strcmp('freq', c.Properties.VariableNames))
    c.fLow  = c.freq(:,1);
    c.fHigh = c.freq(:,2);
end

button = 1;
while button ~= 27
    hTile = tiledlayout(sp.nrow,sp.ncol);

    nShown = min(sp.nrow*sp.ncol, nDet-k+1);
    for i = 1:nShown
        nexttile
        ix = k+i-1;

        switch sp.windowType
            case 'fixed'
                midpoint  = mean([c.t0(ix), c.tEnd(ix)]);
                startTime = midpoint - 0.5*sp.windowDur/86400;
                endTime   = midpoint + 0.5*sp.windowDur/86400;
            otherwise % 'prePost'
                startTime = c.t0(ix)   - sp.pre/86400;
                endTime   = c.tEnd(ix) + sp.post/86400;
        end

        wav = getAudioFromFiles(wavInfo, startTime, endTime, channel=1, newRate=sp.sampleRate);
        [~,f,t,p] = spectrogram(wav, sp.nfft, sp.noverlap, sp.nfft, sp.sampleRate, 'yAxis');

        pdB   = 10*log10(p);
        fIx   = double(f > sp.lowFreq & f < sp.highFreq);
        fIx(fIx==0) = nan;
        pIx   = repmat(fIx,1,length(t));
        h = pcolor(t,f,pdB.*pIx);
        set(h,'lineStyle','none');
        clim = caxis;
        caxis([-50 0]+clim(2));

        if hasFreq
            hold on;
            xLeft  = (c.t0(ix)   - startTime) * 86400;
            xWidth = (c.tEnd(ix) - c.t0(ix))  * 86400;
            rectangle('Position', [xLeft, c.fLow(ix), xWidth, c.fHigh(ix)-c.fLow(ix)], ...
                'EdgeColor','c', 'LineStyle','--', 'LineWidth', 1.2);
            hold off;
        end

        title(sprintf('Detection #%g', ix));
        xlabel(datestr(c.t0(ix)));

        infoStr = '';
        if hasSnr; infoStr = [infoStr sprintf('SNR %.1f dB  ', c.snr(ix))]; end
        if hasRL;  infoStr = [infoStr sprintf('RL %.1f dB  ',  c.RL(ix))];  end
        if hasNL;  infoStr = [infoStr sprintf('NL %.1f dB',    c.NL(ix))];  end
        if ~isempty(infoStr)
            text(0.02,0.95,infoStr,'Units','normalized','Color','w', ...
                'VerticalAlignment','top','FontSize',8);
        end
    end
    pageIx = k:(k+nShown-1);

    [~,~,button] = ginput(1);
    set(hTile.Title,'String','Seeking');
    drawnow;
    switch button
        case 27 % Esc -- quit
            set(hTile.Title,'String','Done');
        case 28 % Left -- previous page
            k = k - sp.nrow*sp.ncol;
        case 29 % Right -- next page
            k = k + sp.nrow*sp.ncol;
        case 30 % Up -- jump forward ~5%
            k = k + floor(nDet*0.05);
        case 31 % Down -- jump back ~5%
            k = k - floor(nDet*0.05);
        case 32 % Space -- redraw
        case 105 % i key -- print current page's detections to console
            qcPrintDetections(c(pageIx,:), 'n', numel(pageIx));
        case 115 % s key -- save current page as PNG
            saveAsPng(sprintf('qcPageDetections_%03g', k));
        otherwise
            k = k + sp.nrow*sp.ncol;
    end

    k = rem(k,nDet);
    if k < 1
        k = nDet + k;
    end
end
end
