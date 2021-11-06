function [h] = plot_line_yint(val,varargin)
% [h] = plot_line_yint([axes handle],val,varargin)
%
% plots a line that extends for a long way along the x axis
%  if given a vector, plots one line for each entry
%
% 'for a long way' means 'until 10^12'
%
% 2019-04-08 edit: optionally, the first argument can be an axes where this
% will be plotted
%
% code by ESBM, 2009, 2019

if isa(val,'matlab.graphics.axis.Axes') || isa(val,'matlab.ui.control.UIAxes')
    assert(~isempty(varargin),'plot_line_yint was given an axes but no data to plot on it!');
    ax = val;
    val = varargin{1};
    varargin = varargin(2:end);
else
    ax = gca;
end
% This code is because oddly, UIFigures cannot handle 
% as long of a line as normal axes...
if matlab.ui.internal.isUIFigure(ancestor(ax,'figure'))
    extent = 10^5;
else
    extent = 10^12;
end

val = torowvector(val);
h = plot(ax,extent.*repmat([-1;1],1,numel(val)),repmat(val,2,1),varargin{:},'xliminclude','off');
