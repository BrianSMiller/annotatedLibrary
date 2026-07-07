function qcDetectionRatePlot(detTable, varargin)
% qcDetectionRatePlot  Coarse pattern-plausibility check: detections per
% day, faceted by site (and optionally by detector/tag), to catch bugs
% that show up as implausible temporal patterns rather than per-detection
% errors -- e.g. a detector firing constantly during a known equipment
% gap, or a suspiciously flat rate across a whole deployment.
%
%   qcDetectionRatePlot(detTable)
%   qcDetectionRatePlot(detTable,'siteCol','site','dateCol','datetime', ...
%       'groupCol','Tags','binWidth',1)
%
% Inputs
%   detTable  - table of detections with a datetime column
%
% Name-value args
%   dateCol    - datetime column name (default 'datetime')
%   siteCol    - site column name (default 'site'); pass '' to skip faceting
%   groupCol   - optional secondary grouping e.g. detector/tag; pass '' to skip
%   binWidth   - histogram bin width in days (default 1)
%   flagCol    - optional column marking detection status. Two modes:
%                logical (e.g. 'qcValid') overlays a single red 'x' per
%                day for flagged counts. String/categorical (e.g. the
%                'qcFlag' column from qcValidateDetections, which is ''
%                for valid rows and a category name otherwise) overlays a
%                distinct colour/marker per category on a secondary
%                right-hand axis -- lets you see whether a specific
%                problem clusters on particular days (a stronger signal
%                of something specific happening then) versus scattering
%                evenly across the whole deployment (more likely a
%                background rate). Passing 'qcFlag' is usually more
%                informative than 'qcValid' for this reason.
%
% B. Miller, AAD, 2026 -- pulled out of exploreDetections.m

p = inputParser;
addParameter(p,'dateCol','datetime');
addParameter(p,'siteCol','site');
addParameter(p,'groupCol','');
addParameter(p,'flagCol','');
addParameter(p,'binWidth',1);
parse(p,varargin{:});
opt = p.Results;

dt = detTable.(opt.dateCol);
edges = dateshift(min(dt),'start','day'):opt.binWidth:dateshift(max(dt),'start','day')+1;
dayCenters = edges(1:end-1) + diff(edges)/2;

hasFlagCol = ~isempty(opt.flagCol) && any(strcmp(opt.flagCol,detTable.Properties.VariableNames));

if hasFlagCol
    flagVals = detTable.(opt.flagCol);
    if islogical(flagVals)
        flagCategories = {'Flagged'};
    else
        flagVals = string(flagVals);
        flagCategories = cellstr(unique(flagVals(flagVals ~= "")));
    end
    markers = {'x','+','*','.','x','+','*','.'};
    flagColors = autumn(max(numel(flagCategories),1));
end

if ~isempty(opt.siteCol) && any(strcmp(opt.siteCol,detTable.Properties.VariableNames))
    sites = unique(detTable.(opt.siteCol));
else
    sites = {'all'};
end

figure;
tl = tiledlayout(length(sites),1);
ax = gobjects(length(sites),1);
for i = 1:length(sites)
    ax(i) = nexttile;
    if isequal(sites,{'all'})
        rows = true(height(detTable),1);
    else
        rows = strcmp(cellstr(string(detTable.(opt.siteCol))),string(sites(i)));
    end
    if ~isempty(opt.groupCol) && any(strcmp(opt.groupCol,detTable.Properties.VariableNames))
        groups = unique(detTable.(opt.groupCol)(rows));
        hold on;
        for g = 1:length(groups)
            gRows = rows & strcmp(cellstr(string(detTable.(opt.groupCol))),string(groups(g)));
            histogram(dt(gRows),edges,'DisplayName',string(groups(g)));
        end
        legend;
        hold off;
    else
        histogram(dt(rows),edges);
    end
    ylabel('Detections/day');
    title(string(sites(i)),'Interpreter','none');
    grid on;

    if hasFlagCol
        yyaxis right
        hold on;
        for c = 1:numel(flagCategories)
            if islogical(flagVals)
                catMask = rows & ~flagVals;
            else
                catMask = rows & (flagVals == flagCategories{c});
            end
            catCounts = histcounts(dt(catMask), edges);
            hasAny = catCounts > 0;
            if any(hasAny)
                plot(dayCenters(hasAny), catCounts(hasAny), ...
                    'Marker', markers{mod(c-1,numel(markers))+1}, ...
                    'LineStyle', 'none', 'Color', flagColors(c,:), ...
                    'MarkerSize', 7, 'LineWidth', 1.3, ...
                    'DisplayName', flagCategories{c});
            end
        end
        hold off;
        ylabel('Flagged/day');
        if numel(flagCategories) > 1
            legend('Location','eastoutside');
        end
        yyaxis left
    end
end
linkaxes(ax,'x');
tl.XLabel.String = 'Date';
end
