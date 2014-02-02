% draw_line(x, y, k, xmin, xmax, style)
%
% Draws a line from (x, y) with slope k between
% xmin and xmax.
function draw_line(x, y, k, xmin, xmax, style)
	ymin = y + (xmin - x) * k;
	ymax = y + (xmax - x) * k;
	plot([xmin xmax], [ymin ymax], style)
end