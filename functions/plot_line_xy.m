function [h] = plot_line_xy(x,y,varargin)
% [h] = plot_line_xy([axes handle],x,y,varargin)
%
% plots a line that extends for a long way
% and passes through the two given (x,y) points
%
% 'for a long way' means 'until 10^7 from the first point'
%
% 2016-06-10 edit: changed scaling factor from 10^12 to 10^5, because
%  the new version of matlab doesn't seem to handle plotting
%  big numbers with high precision properly!
%
% 2019-04-08 edit: optionally, the first argument can be an axes where this
% will be plotted
%
% code by ESBM, 2007, 2016, 2019

if (isa(x,'matlab.graphics.axis.Axes') || isa(x,'matlab.ui.control.UIAxes'))
    assert(~isempty(varargin),'plot_line_xy was given an axes but no data to plot on it!');
    ax = x;
    x = y;
    y = varargin{1};
    varargin = varargin(2:end);
else
    ax = gca;
end
% This code is because oddly, UIFigures cannot handle 
% as long of a line as normal axes...
if matlab.ui.internal.isUIFigure(ancestor(ax,'figure'))
    extent = 10^4;
else
    extent = 10^5;
end

if ~isvector(x) || ~isvector(y) || length(x) ~= 2 || length(y) ~= 2
    error('inputs must be 2-length vectors');
end

normvect = normalize([diff(x) diff(y)],2);

h = plot(ax, ...
    x(1) + extent*normvect(1)*[-1 1], ...
    y(1) + extent*normvect(2)*[-1 1], ...
    varargin{:},'xliminclude','off','yliminclude','off');
