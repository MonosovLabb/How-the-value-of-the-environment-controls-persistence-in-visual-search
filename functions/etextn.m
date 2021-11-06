function h = etextn(alignsettings,varargin)
% h = etextn([axes,]alignsettings,varargin)
% call 'text' with Ethan's usual settings
%  using horizontalalign = ha, verticalalign = va
%
% 'etextn' is the same as 'etext', except for two things:
%
%  1. It applies the setting: 'units','normalized' 
%
%  2. You can choose not to specify the x and y coordinates.
%     If your first argument after the alignment settings, and after any
%     leading entries in varargin that specify the plotting axes,
%     is not numeric then, it assumes that you want x and y
%     to be chosen automatically. Then, sets them to:
%
%      x: 0.5, 0.02, 0.98 (for horizontal: c, l, r)
%      y: 0.5, 0.98, 0.02, 0.98, 0.02 (for vertical: m, t, b, c, l)
%
% alignsettings:
%  horizontal: (c)enter, (l)eft, (r)ight
%  vertical:   (m)iddle, (t)op, (b)ottom, (c)ap, base(l)ine
%  (default = 'cm')
%
% Code by ESBM, 2010, 2020

varargin = [varargin 'units','normalized'];

% get axis to be used to plot the data: either a user-specified axis,
% or default to 
if numel(varargin) >= 1 && (isa(alignsettings,'matlab.graphics.axis.Axes') || isa(alignsettings,'matlab.ui.control.UIAxes'))
    ax = alignsettings;
    alignsettings = varargin{1};
    varargin = varargin(2:end);
else
    ax = gca;
end

% assume that user specified the x,y coordinates if the first argument is
% numeric, or the first 
user_specified_xy_coords = numel(varargin) >= 1 && isnumeric(varargin{1});


if ~user_specified_xy_coords
    % set x, y automatially!
    
    if isempty(alignsettings)
        alignsettings = 'cm';
    end

    switch alignsettings(1)
        case 'c'
            x = 0.5;
        case 'l'
            x = 0.01;
        case 'r'
            x = 0.99;
        otherwise
            error('invalid horizontal alignment setting');
    end
    switch alignsettings(2)
        case 'm'
            y = 0.5;
        case 't'
            y = 0.99;
        case 'b'
            y = 0.01;
        case 'c'
            y = 0.99;
        case 'l'
            y = 0.01;
        otherwise
            error('invalid vertical alignment setting');
    end
    
    h = etext(ax,alignsettings,x,y,varargin{:});
else
    % use user specified x,y
    h = etext(ax,alignsettings,varargin{:});
end
