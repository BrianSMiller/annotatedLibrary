function qcSnrDiagnosticPlots(detTable, varargin)
% qcSnrDiagnosticPlots  SNR/NL/RL diagnostic figure: monthly time series
% (boxplot) for each of SNR, NL, RL, paired with a histogram rotated to
% share the same y-axis scale -- lets you read distribution shape and
% temporal pattern side by side rather than in two separate figures.
%
%   qcSnrDiagnosticPlots(detTable)
%   qcSnrDiagnosticPlots(detTable, 'dateCol','datetime', 'widthRatio', 3)
%
% Inputs
%   detTable  - detection table with snr, NL, RL (from qcAttachSnr) and
%               a datetime column
%
% Name-value args
%   dateCol     - datetime column name (default 'datetime')
%   widthRatio  - relative width of time-series panel vs histogram panel
%                 (default 4, i.e. 4:1)
%   validOnly   - if true (default), restricts to qcValid rows when that
%                 column exists -- SNR isn't meaningful for flagged rows
%
% Boxplot ('PlotStyle','compact') stands in for a violin plot -- MATLAB
% has no built-in violin. If you have one on your path (gramm, or a
% FileExchange violin.m), the month-grouped data here is ready to hand
% straight to it instead.
%
% B. Miller, AAD, 2026

p = inputParser;
addParameter(p,'dateCol','datetime');
addParameter(p,'widthRatio',4);
addParameter(p,'validOnly',true);
parse(p,varargin{:});
opt = p.Results;

if opt.validOnly && any(strcmp('qcValid',detTable.Properties.VariableNames))
    detTable = detTable(detTable.qcValid, :);
end
vars   = {'snr','NL','RL'};
titles = {'SNR (dB)', 'NL (dB re 1 \muPa)', 'RL (dB re 1 \muPa)'};

warning('off', 'stats:boxplot:DeprecatedTiledChartLayout');
warning('off', 'stats:boxplot:BadObjectType')

nCols = opt.widthRatio + 1;
figure;
tl = tiledlayout(3, nCols, 'TileSpacing','compact');

dt = detTable.(opt.dateCol);

for i = 1:numel(vars)
    varName = vars{i};
    if ~any(strcmp(varName, detTable.Properties.VariableNames))
        continue;
    end
    x = detTable.(varName);
    finiteMask = isfinite(x);
    if ~any(finiteMask)
        nexttile([1 nCols]);
        title(sprintf('%s (no finite values)', titles{i}));
        continue;
    end

    monthsDT    = dateshift(dt(finiteMask), 'start', 'month');
    monthLabels = cellstr(datestr(monthsDT, 'yyyy-mm'));

    % Time series panel (spans widthRatio of nCols columns)
    ax1 = nexttile([1 opt.widthRatio]);

    boxplot(x(finiteMask), monthLabels, 'PlotStyle','compact');
    ylabel(titles{i});
    title(sprintf('%s -- monthly distribution', titles{i}));
    xtickangle(45);
    grid on;
    yl = ylim(ax1);

    % Histogram panel (1 column), sharing y-limits with the time series
    ax2 = nexttile([1 1]);
    histogram(x(finiteMask), 40, 'Orientation','horizontal', ...
        'FaceColor',[0.2 0.4 0.8], 'FaceAlpha',0.7);
    ylim(ax2, yl);
    xlabel('Count');
    set(ax2, 'YTickLabel', []);
    grid on;

    
end

tl.Title.String = 'SNR / NL / RL diagnostics';
tl.Title.Interpreter='none';

% Add after the existing vars loop in qcSnrDiagnosticPlots, before the
% tl.Title line. Requires RL and (fLow,fHigh) present.

hasBandCols = all(ismember({'RL','fLow','fHigh'}, detTable.Properties.VariableNames));
if hasBandCols
    bw = detTable.fHigh - detTable.fLow;
    rl = detTable.RL;
    good = isfinite(bw) & isfinite(rl) & bw > 0;

    if any(good)
        % Full-width row for the scatter
        figure;
        scatter(bw(good), rl(good), 8, 'filled', ...
            'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', [0.2 0.4 0.8]);
        xlabel('Bandwidth (Hz)');
        ylabel('RL band level (dB re 1 \muPa)');
        title('RL vs bandwidth (correlation \rightarrow low levels are narrow, not quiet)');
        grid on;

        % Overlay the correlation and a robust fit line so the relationship
        % is quantified, not just eyeballed
        r = corr(bw(good), rl(good), 'type', 'Spearman');
        hold on;
        pf = polyfit(bw(good), rl(good), 1);
        xl = xlim;
        plot(xl, polyval(pf, xl), 'r-', 'LineWidth', 1.2);
        text(0.02, 0.95, sprintf('Spearman \\rho = %.2f', r), ...
            'Units','normalized', 'VerticalAlignment','top', ...
            'BackgroundColor','w', 'FontSize', 9);
        hold off;
    end
end


end
