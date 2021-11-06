function [h] = liney(val,varargin)
% [h] = liney(val,varargin)
%
% plots a line that extends for a long way along the x axis
%  if given a vector, plots one line for each entry
%
% 'for a long way' means 'until 10^12'
% 
% code by ESBM, 2009

h = plot_line_yint(val,varargin{:});
