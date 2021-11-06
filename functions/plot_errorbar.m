function [h] = plot_errorbar(x,y,lo,hi,varargin)
% [h] = plot_errorbar(x,y,lo,hi,varargin)
% draw vertical lines through each point to make error bars
%
% if varargin gets two cell arrays as arguments, uses the first
%  as the style for the main plot line, and the second as
%  the style for the error bars.
% otherwise, variable arguments are used for properties of both.
%
% if the first variable argument for the main plot line is 'bar',
%  then is plotted as a separate bar for each data point, and
%  then second variable argument is taken as bar face color.
%
% code by ESBM, 2011, 2017

if ~isvector(x) || ~isvector(y)
    error('inputs must be vectors');
end;

if nargin < 3
    error('too few args');
elseif nargin == 3
    hi = y + lo;
    lo = y - lo;
elseif nargin >= 4
    if ischar(hi) || iscell(hi)
        varargin = [{hi} varargin];
        hi = y + lo;
        lo = y - lo;
    end;
end;

x = x(:)';
y = y(:)';
lo = lo(:)';
hi = hi(:)';

use_separate_styles = length(varargin) == 2 && iscell(varargin{1}) && iscell(varargin{2});
if use_separate_styles
    arg1 = varargin{1};
    arg2 = varargin{2};
else
    arg1 = varargin;
    arg2 = varargin;
end;

was_hold = ishold;
hold on;

% handle 'bar' style plots
if numel(arg2) >= 2 && isequal(arg2{1},'bar')
    % ignore 'bar' arguments for errorbars
    arg2 = arg2(3:end);
end;
if numel(arg1) >=2 && isequal(arg1{1},'bar')
    bar_facecolor = arg1{2};
    bar_args = arg1(3:end);
    h1 = plot_bar(x,y,bar_facecolor,bar_args{:});
    h2 = plot([x ; x],[lo ; hi],arg2{:});
else
    % normal plotting style
    h2 = plot([x ; x],[lo ; hi],arg2{:});
    h1 = plot(x,y,arg1{:});
end;


h = [h1 ; h2];

if ~was_hold
    hold off;
end;

