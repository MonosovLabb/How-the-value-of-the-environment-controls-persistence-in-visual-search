function [h] = figuren(varargin)
% [h] = figuren(varargin)
% same as 'figure', but applies my favorite settings using 'nsubplot'!
%
% Code by ESBM, 2012

h = figure(varargin{:});set(gcf,'color',[1 1 1])
nsubplot;