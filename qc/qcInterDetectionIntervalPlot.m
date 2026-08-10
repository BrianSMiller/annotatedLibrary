function stats = qcInterDetectionIntervalPlot(detTable, expectedRange_s, varargin)
% qcInterDetectionIntervalPlot  Inter-detection-interval (IDI) diagnostic.
%
% Species-agnostic single-observer QC check: compares the distribution of
% intervals between successive detections against an expected biological
% range for the target call type. Same underlying idea used for sperm
% whale click-train classification (spermWhaleIndex.m/plotIdiGram.m) and
% ABW song IDI in Common Ground -- generalised here with no species-
% specific limits baked in, run once per detector per site before
% matching.
%
%   qcInterDetectionIntervalPlot(detTable, expectedRange_s)
%   qcInterDetectionIntervalPlot(detTable, [8 20], 'sliceDuration_s', 1800)
%   stats = qcInterDetectionIntervalPlot(...)
%
% Inputs
%   detTable         detection table with a datetime column (added by
%                     qcDetector.m before this runs)
%   expectedRange_s  [low high] expected IDI range in seconds for the
%                     target call type, e.g.:
%                       sperm whale clicks   [0.5  2]
%                       fin whale 20 Hz      [8   20]
%                       ABW song             [60  70]
%                     No default -- must be supplied per call type.
%                     Must satisfy high < sliceDuration_s.
%
% Name-value args
%   dateCol          - datetime column name (default 'datetime')
%   validCol         - logical column restricting to valid rows when it
%                      exists (default 'qcValid'); pass '' to use all rows
%   nIdiBins         - log-spaced IDI histogram bins (default 60)
%   sliceDuration_s  - frame/slice duration in seconds (default 3600).
%                      Serves two roles: (1) the time-resolved panel's
%                      bin width, and (2) the forced upper edge of the
%                      IDI histogram bins -- any IDI >= sliceDuration_s
%                      is clipped into the top bin for BINNING purposes
%                      only. Deployment-scale gaps (equipment outages,
%                      duty-cycle silence) are already covered by
%                      qcDetectionRatePlot; letting them stretch this
%                      histogram's log axis out to days/weeks would
%                      compress the biologically-relevant range this
%                      check actually exists to look at. stats.idi_s
%                      always holds the true, unclipped intervals, so
%                      pctAboveRange etc. are unaffected by the clip.
%
% Detections are assumed unable to occur in gaps between recording files
% -- if they regularly do, that's a detector/config bug to fix at the
% source (spectrogram settings drifted from the detector's own, a
% timestamp format mismatch, or a corrupted file handled differently by
% the detector than by soundFolder), not something to paper over with
% per-file grouping here. IDI is a single sorted sequence across the
% whole input, no per-file breaks.
%
% Output
%   stats  struct with idi_s (true, unclipped), pctInRange/
%          pctBelowRange/pctAboveRange, medianIdi_s (all computed on the
%          unclipped intervals), plus idiBinEdges/idiGram (clipped, as
%          plotted).
%
% Two panels, same time-series+histogram grammar as qcSnrDiagnosticPlots:
%   Left   - IDI over time, log-spaced y-axis, expected range marked.
%            Adapted from spermWhaleIndex.m/plotIdiGram.m's idiGram --
%            lets a persistent problem (e.g. a systematic multipath
%            artefact) be told apart from an episodic one (e.g. a chorus
%            period, or a genuine gap in calling).
%   Right  - aggregate IDI histogram (log x-axis), expected range shaded,
%            with 2x/3x-low-bound reference lines -- a peak there usually
%            means every-other/every-third pulse is being missed, not
%            that the call type has a different natural spacing.
%
% B. Miller, AAD, 2026

p = inputParser;
addParameter(p,'dateCol','datetime');
addParameter(p,'validCol','qcValid');
addParameter(p,'nIdiBins',60);
addParameter(p,'sliceDuration_s',3600);
parse(p,varargin{:});
opt = p.Results;

if numel(expectedRange_s)~=2 || expectedRange_s(1)>=expectedRange_s(2)
    error('qcInterDetectionIntervalPlot:badRange', ...
        'expectedRange_s must be [low high] with low < high.');
end
if expectedRange_s(2) >= opt.sliceDuration_s
    error('qcInterDetectionIntervalPlot:rangeExceedsSlice', ...
        ['expectedRange_s upper bound (%.4g s) must be less than ' ...
         'sliceDuration_s (%.4g s).'], expectedRange_s(2), opt.sliceDuration_s);
end

if ~isempty(opt.validCol) && any(strcmp(opt.validCol,detTable.Properties.VariableNames))
    detTable = detTable(detTable.(opt.validCol), :);
end

dt = sort(detTable.(opt.dateCol));
if numel(dt) < 2
    warning('qcInterDetectionIntervalPlot:tooFewDetections', ...
        'Fewer than 2 valid detections -- nothing to compute.');
    stats = struct('idi_s',[],'pctInRange',NaN,'pctBelowRange',NaN, ...
        'pctAboveRange',NaN,'medianIdi_s',NaN);
    return;
end

idi_s   = seconds(diff(dt));          % true, unclipped intervals
idiTime = dt(2:end);                  % IDI attributed to the later detection in each pair

stats.idi_s         = idi_s;
stats.pctInRange    = 100*mean(idi_s>=expectedRange_s(1) & idi_s<=expectedRange_s(2));
stats.pctBelowRange = 100*mean(idi_s<expectedRange_s(1));
stats.pctAboveRange = 100*mean(idi_s>expectedRange_s(2));
stats.medianIdi_s   = median(idi_s);

% Clip only for binning/plotting -- stats above stay on the true values.
idiClipped = min(idi_s, opt.sliceDuration_s);

%% Time-resolved idiGram: log-spaced IDI bins x sliceDuration_s time bins
minEdge = max(min(idi_s),1e-3);
if minEdge >= opt.sliceDuration_s
    minEdge = opt.sliceDuration_s/1e3; % degenerate-data fallback
end
idiBinEdges   = logspace(log10(minEdge), log10(opt.sliceDuration_s), opt.nIdiBins+1);
idiBinCentres = sqrt(idiBinEdges(1:end-1).*idiBinEdges(2:end)); % geometric centre, log bins

timeBinWidth = seconds(opt.sliceDuration_s);
timeStart = dateshift(min(idiTime),'start','hour');
timeEnd   = dateshift(max(idiTime),'start','hour') + timeBinWidth;
timeEdges = timeStart:timeBinWidth:timeEnd;
nTimeBins = numel(timeEdges)-1;

idiGram = zeros(opt.nIdiBins, nTimeBins);
for i = 1:nTimeBins
    inBin = idiTime>=timeEdges(i) & idiTime<timeEdges(i+1);
    if any(inBin)
        idiGram(:,i) = histcounts(idiClipped(inBin), idiBinEdges);
    end
end

stats.timeEdges     = timeEdges;
stats.idiBinEdges   = idiBinEdges;
stats.idiBinCentres = idiBinCentres;
stats.idiGram       = idiGram;

%% Plot
figure;
tl = tiledlayout(1,4,'TileSpacing','compact');

hIdiGram = nexttile([1 3]);
tbNum = datenum(timeEdges(1:end-1)); %#ok<DATNM> -- imagesc needs numeric XData
hpcol = pcolor(tbNum, idiBinCentres, idiGram);
set(hpcol,'LineStyle','none')
set(gca,'YDir','normal','YScale','log');
datetick('x','keeplimits');
yt = get(gca,'YTick');
% set(gca,'YTickLabel', compose('%.2g', 10.^yt));
colormap(flipud(gray));
colorbar('north');
ylabel('Inter-detection interval (s)');
xlabel('Time');
hold on;
yline(expectedRange_s(1),'r--','LineWidth',1);
yline(expectedRange_s(2),'r--','LineWidth',1);
grid on;
hold off;
title(sprintf('IDI over time (%g s slices), expected range marked', opt.sliceDuration_s));

hIdiCount = nexttile;
histogram(idiClipped, idiBinEdges,'orientation','horizontal');
set(gca,'YScale','log');
hold on;
yl = ylim;
linkaxes([hIdiCount, hIdiGram],'y')
ylim([minEdge,opt.sliceDuration_s]);
patch([yl(1) yl(1) yl(2) yl(2)],...
    [expectedRange_s(1) expectedRange_s(2) expectedRange_s(2) expectedRange_s(1)],...
    [0.9 0.9 0.9], 'FaceAlpha',0.4, 'EdgeColor','r');
yline(stats.medianIdi_s,'m-','median','LineWidth',1.5,'LabelVerticalAlignment','middle');
yline(2*expectedRange_s(1),'k--','2x low bound','LabelVerticalAlignment','middle');
yline(3*expectedRange_s(1),'k:','3x low bound','LabelVerticalAlignment','middle');
uistack(findobj(gca,'Type','histogram'),'top');
xlabel('Inter-detection interval (s, log scale)');
ylabel('Count');
text(max(xlim),max(ylim),sprintf(['%.0f%% in [%.1f %.1f]s, \n',...
    '%.0f%% below, %.0f%% above \n', ...
    '(>= %g s clipped to top bin)'], ...
    stats.pctInRange, expectedRange_s(1), expectedRange_s(2), ...
    stats.pctBelowRange, stats.pctAboveRange, opt.sliceDuration_s), ...
    'VerticalAlignment','top','HorizontalAlignment','right');
hold off;

tl.Title.String = 'Inter-detection interval diagnostic';
tl.Title.Interpreter = 'none';

end
