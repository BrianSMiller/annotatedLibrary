function [tSpan, fSpan] = drawBox(ax, firstPoint, color)
% DRAWBOX  Draw a box on a spectrogram with two clicks at opposite corners
%
%   [tSpan, fSpan] = drawBox(ax) waits for two clicks in axes ax and returns
%   the box as sorted [start end] and [low high] in the axes' units (e.g.
%   seconds from the start of the window, and Hz). The box is drawn green.
%
%   [tSpan, fSpan] = drawBox(ax, firstPoint) takes the first corner as
%   already clicked ([x y], e.g. from ginput) and waits for the second.
%
%   [tSpan, fSpan] = drawBox(ax, firstPoint, color) draws in that colour.
%
%   Any mouse button places a corner. Esc (or any key) before the second
%   corner cancels: tSpan and fSpan are returned empty.
%
%   Used by judgeDetections (the R key, with the 'rebox' option) and by
%   reboxTagCalls (awr-covariates).
%
%   See also JUDGEDETECTIONS, GINPUT

arguments
    ax = gca
    firstPoint (1,:) double = []
    color = 'g'
end

tSpan = [];  fSpan = [];
axes(ax);
if isempty(firstPoint)
    [x1, y1, b1] = ginput(1);
    if isempty(b1) || ~ismember(b1, 1:3), return; end
else
    x1 = firstPoint(1);  y1 = firstPoint(2);
end
wasHeld = ishold(ax);
hold(ax, 'on');
restoreHold = onCleanup(@() restore(ax, wasHeld));   % leave hold as it was
hMark = plot(ax, x1, y1, '+', 'Color', color, 'MarkerSize', 12, 'LineWidth', 1.5);
drawnow;
[x2, y2, b2] = ginput(1);
delete(hMark);
if isempty(b2) || ~ismember(b2, 1:3), return; end

tSpan = sort([x1 x2]);
fSpan = sort([y1 y2]);
plot(ax, tSpan([1 1 2 2 1]), fSpan([1 2 2 1 1]), 'Color', color, 'LineWidth', 1.5);
drawnow;
end

function restore(ax, wasHeld)
if ~wasHeld && isvalid(ax), hold(ax, 'off'); end
end
