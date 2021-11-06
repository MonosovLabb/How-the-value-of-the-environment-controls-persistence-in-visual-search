function [h] = linexy(varargin)
% [h] = linexy(...)
%
%  plots a line that extends for a long way
%  and passes through (0,0) and (1,1)
% 
% [h] = linexy(x,y,...)
%
%  plots a line that extends for a long way
%  and passes through the two given (x,y) points
%
% 'for a long way' means 'until 10^12 from the first point'
% 
% code by ESBM, 2009

if numel(varargin) >=1 && (isa(varargin{1},'matlab.graphics.axis.Axes') || isa(varargin{1},'matlab.ui.control.UIAxes'))
    ax = varargin{1};
    varargin = varargin(2:end);
else
    ax = gca;
end

if numel(varargin) < 2 || ~isnumeric(varargin{1})
    h = plot_line_xy(ax,[0 1],[0 1],varargin{:});
else
    h = plot_line_xy(ax,varargin{:});
end
