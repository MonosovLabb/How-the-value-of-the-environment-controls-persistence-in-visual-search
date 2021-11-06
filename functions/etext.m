function h = etext(alignsettings,varargin)
% h = etext([axes],alignsettings,varargin)
% call 'text' with Ethan's usual settings
%  using horizontalalign = ha, verticalalign = va
%
% alignsettings:
%  horizontal: (c)enter, (l)eft, (r)ight
%  vertical:   (m)iddle, (t)op, (b)ottom, (c)ap, base(l)ine
%  (default = 'cm')
%
% Code by ESBM, 2010, 2020

if numel(varargin) >= 1 && (isa(alignsettings,'matlab.graphics.axis.Axes') || isa(alignsettings,'matlab.ui.control.UIAxes'))
    ax = alignsettings;
    alignsettings = varargin{1};
    varargin = varargin(2:end);
else
    ax = gca;
end

if isempty(alignsettings)
    alignsettings = 'cm';
end

ha = 'center';
va = 'middle';
switch alignsettings(1)
    case 'c'
        ha = 'center';
    case 'l'
        ha = 'left';
    case 'r'
        ha = 'right';
    otherwise
        error('invalid horizontal alignment setting');
end
switch alignsettings(2)
    case 'm'
        va = 'middle';
    case 't'
        va = 'top';
    case 'b'
        va = 'bottom';
    case 'c'
        va = 'cap';
    case 'l'
        va = 'baseline';
    otherwise
        error('invalid vertical alignment setting');
end

spec_units = false;
for v = 1:numel(varargin)
    if ischar(varargin{v}) && ~isempty(strmatch('unit',varargin{v}))
        spec_units = true;
        break;
    end
end

if ~spec_units
    varargin = [varargin {'units','data'}];
end
varargin = [varargin {'horizontalalignment',ha,'verticalalignment',va}];

h = text(ax,varargin{:});
