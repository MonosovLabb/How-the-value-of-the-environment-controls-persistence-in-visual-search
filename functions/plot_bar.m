function [h] = plot_bar(x,n,color,varargin)
% [h] = plot_bar(x,n,color,...)
%
% bar plot
%  x = x coords
%   (can be either the midpoints of the bars, or the edges)
%  n = height of the bars
% ... = arguments passed to plotting function
%  but can start with up to four special arguments:
% 
%   (any scalar) = width of the bins on the plot (1=max width, <1 = leave empty space)
%   'fixedwidth',(any scalar) = force all bins to be a specified width (NOT a relative width, as above)
%   'vert' = plot bars vertically instead of horizontally
%   'line' = plot using lines instead of bars
%   'normalize' = normalize n values to sum to 1
%   'normalize',N = normalize n values to sum to N
%   'negative' = multiply bar heights by -1 before plotting
%   'outline',outlinecolor = plot with an outline around the histogram using
%       the specified color (or, if "outlinecolor" is a cell array, use the
%       arguments given in the cell array)
%   'barbottom',B = bars start at y=B (default = bars start at y=0)
%     NOTE: deprecated: this was implemented by simply adding a vertical
%     offset to the entire bar, instead of just changing its start
%     position! Use the below argument instead:
%   'barstart',B = bars start at y=B (default = 0)
%   'axes',A = plot to axes A (default = the current axis)
%   'mergebars' = merge adjacent patch bars into single 'patch' objects
%     (useful for e.g. minimizing file size in .eps file)
%
%  if first argument is an axis, plots to that axis...
%
% example:
%  plot_bar(xvals,counts,'r',.8,'vert','normalize','EdgeColor','b')
%
% warning: do not use directly on output of 'histc'
%  but instead, remove the last 'count' from the data.
%  because histc includes a special 'count' for data
%  exactly equal to the last edge value, which messes
%  up this type of plot!
%
% code by ESBM, 2012


% if no data, don't plot it!
origx = x;
orign = n;

fixed_bar_width = nan;
drawn_bin_width = 1;
plot_vertically = 0;
plot_with_lines = 0;
plot_negative_bars = false;
outlinecolor = nan;
normalize_values = 'no';
barbottom = nan;
barstart = nan;
plotaxes = gca;
mergebars = false;

if nargin < 3 || isempty(color)
    color = [.7 .7 .7];
end;
if any(ischar(color))
    switch color
        case {'vert','line','normalize','negative','barbottom','barstart','axes','outline'}
            varargin = [{color} varargin];
            color = [.7 .7 .7];
    end;
end;

% read arguments
v = 1;
while true
    if v > length(varargin)
        break;
    end;
    
    if isscalar(varargin{v}) && ~ischar(varargin{v})
        drawn_bin_width = varargin{v};
        v = v + 1;
    elseif strcmp(varargin{v},'fixedwidth')
        if length(varargin) >= v+1
            fixed_bar_width = varargin{v+1};
            v = v + 2;
        else
            error('expected a width after ''fixedwidth''');
        end;
    elseif strcmp(varargin{v},'vert')
        plot_vertically = true;
        v = v + 1;
    elseif strcmp(varargin{v},'line')
        plot_with_lines = true;
        v = v + 1;
    elseif strcmp(varargin{v},'normalize')
        if length(varargin) >= v+1 && isscalar(varargin{v+1})
            normalize_values = varargin{v+1};
            v = v + 2;
        else
            normalize_values = 1;
            v = v + 1;
        end;
    elseif strcmp(varargin{v},'negative')
        plot_negative_bars = true;
        v = v + 1;
    elseif strcmp(varargin{v},'outline')
        if length(varargin) >= v+1
            outlinecolor = varargin{v+1};
            v = v+2;
        else
            error('expected an outline color after ''outline''');
        end;
    elseif strcmp(varargin{v},'barbottom')
        if length(varargin) >= v+1 && isscalar(varargin{v+1})
            barbottom = varargin{v+1};
            v = v+2;
        else
            error('expected value ''B'' after input ''barbottom''');
        end;
    elseif strcmp(varargin{v},'barstart')
        if length(varargin) >= v+1 && isscalar(varargin{v+1})
            barstart = varargin{v+1};
            v = v+2;
        else
            error('expected value ''B'' after input ''barstart''');
        end;
    elseif strcmp(varargin{v},'axes')
        if length(varargin) >= v+1 && isscalar(varargin{v+1})
            plotaxes = varargin{v+1};
            assert(ishandle(plotaxes) && strcmp('axes',get(plotaxes,'type')));
            v = v + 2;
        end;
    elseif strcmp(varargin{v},'mergebars')
        mergebars = true;
        v = v +1;
    else
        break;
    end;
end;
varargin = {varargin{v:end}};

if isnan(barbottom)
    if isnan(barstart)
        barstart = 0;
    end
else
    if isnan(barstart)
        disp(' plot_bar: warning: "barbottom" argument is deprecated, use "barstart"');
    else
        disp(' plot_bar: warning: both "barbottom" and "barstart" arguments were specified as non-NaN, will only use "barstart"');
        barbottom = nan;
    end
end
assert( isnan(barstart) ~= isnan(barbottom),'expected "barstart" or "barbottom" to be non-NaN and the other to be NaN; only one should be specified');


%normalize
if isscalar(normalize_values)
    n = normalize_values.*normalize(n,1);
end;

% figure out where to place bars on the x axis
if length(x) == length(n)+1
    binstart = x(1:(end-1));
    binmid = (x(1:(end-1)) + x(2:end)) ./ 2;
    binend = x(2:end);
elseif length(x) == length(n)
    if length(x) == 1
        binmid = x;
        binstart = binmid-.5;
        binend = binmid+.5;
    else
        binmid = x;
        binstart = (binmid + prev(binmid)) ./ 2;
        binstart(1) = binmid(1) - (binmid(2) - binmid(1))./2;
        binend = (binmid + next(binmid)) ./ 2;
        binend(end) = binmid(end) + (binmid(end) - binmid(end-1))./2;
    end;
else
    error('x must have length equal to length(n) or length(n)+1');
end;
if ~isnan(fixed_bar_width)
    xmin = binmid - fixed_bar_width./2;
    xmax = binmid + fixed_bar_width./2;
else
    xmin = binmid - (binmid - binstart).*drawn_bin_width;
    xmax = binmid + (binend - binmid).*drawn_bin_width;
end;

if ~isnan(barstart)
    ymin = barstart;
else
    % preserve old bad behavior of deprecated "barbottom" 
    % argument, to not break old code
    assert(~isnan(barbottom));
    ymin = 0;
end

% plot
if all(n >= 0)
    patchids = find(n > 0);
else
    patchids = 1:numel(n);
end;

x = [];
y = [];

if plot_with_lines
    for p = 1:length(patchids)
        curx = [];
        cury = [];
        if p == 1 || patchids(p-1) ~= patchids(p)-1
            curx = [curx nan xmin(patchids(p))];
            cury = [cury nan ymin];
        end;
        curx = [curx xmin(patchids(p)) xmax(patchids(p))];
        cury = [cury n(patchids(p)) n(patchids(p))];
        if p == length(patchids) || patchids(p+1) ~= patchids(p)+1
            curx = [curx xmax(patchids(p))];
            cury = [cury ymin];
        end;
        x = [x curx];
        y = [y cury];
    end;
    
    if ~isnan(barstart)
        if plot_negative_bars
            y = -(y-ymin) + ymin;
        end
    else
        assert(~isnan(barbottom));
        % preserve old bad behavior of deprecated "barbottom" 
        % argument, to not break old code
        if plot_negative_bars
            y = barbottom - y;
        else
            y = barbottom + y;
        end;
    end
        
    
    
    if plot_vertically
        h = plot(plotaxes,y,x,'Color',color,varargin{:});
    else
        h = plot(plotaxes,x,y,'Color',color,varargin{:});
    end;
elseif mergebars
    
    okpatch = n > 0;
    cur_patch_start = 1;
    cur_patch_end = 0;
    
    oldaxis = gca;
    axes(plotaxes);
    
    h = [];
    while true
        % find start of next stretch of 'ok' data points
        while cur_patch_start <= numel(okpatch) && ~okpatch(cur_patch_start)
            cur_patch_start = cur_patch_start+1;
        end;
        % no more stretches left
        if cur_patch_start > numel(okpatch)
            break;
        end;
        
        % find end of stretch of 'ok' data points
        cur_patch_end = cur_patch_start;
        while cur_patch_end <= numel(okpatch) && okpatch(cur_patch_end)
            cur_patch_end = cur_patch_end+1;
        end;
        cur_patch_end = cur_patch_end - 1;
        
        x = [];
        y = [];
        for i = cur_patch_start:cur_patch_end
            x = [x binstart(i) binend(i)];
            y = [y n(i) n(i)];
        end;
        if plot_negative_bars
            y = -y;
        end;
        x = [binstart(cur_patch_start) x binend(cur_patch_end) binstart(cur_patch_start)];
        if ~isnan(barstart)
            y = [ymin y ymin ymin];
        else
            % preserve old behavior of deprecated "barbottom" argument, 
            % (which seems to have been correct in this case?)
            % to not break old code
            assert(~isnan(barbottom))
            y = [barbottom y barbottom barbottom];
        end
        
        if plot_vertically
            h_tmp = patch(y,x,color,varargin{:});
        else
            h_tmp = patch(x,y,color,varargin{:});
        end;
        h = [h; h_tmp(:)];
        
        cur_patch_start = cur_patch_end+1;
    end;
    
    axes(oldaxis);
    
else
    curxmin = torowvector(xmin(patchids));
    curxmax = torowvector(xmax(patchids));
    curn = torowvector(n(patchids));
    x = [curxmin ; curxmin ; curxmax ; curxmax];
    y = [ymin*ones(size(curn)) ; curn ; curn; ymin*ones(size(curn))];

    if ~isnan(barstart)
        if plot_negative_bars
            y = -(y-ymin) + ymin;
        end
    else
        assert(~isnan(barbottom));
        % preserve old bad behavior of deprecated "barbottom" 
        % argument, to not break old code
        if plot_negative_bars
            y = barbottom - y;
        else
            y = barbottom + y;
        end;
    end

    oldaxis = gca;
    if oldaxis ~= plotaxes
        axes(plotaxes);
    end
    if plot_vertically
        h = patch(y,x,color,varargin{:});
    else
        h = patch(x,y,color,varargin{:});
    end;
    if oldaxis ~= plotaxes
        axes(oldaxis);
    end
end;

if iscell(outlinecolor) || any(~isnan(outlinecolor))
    extra_args = {};
    if plot_negative_bars
        extra_args{end+1} = 'negative';
    end;
    if plot_vertically
        extra_args{end+1} = 'vert';
    end;
    if isscalar(normalize_values)
        extra_args{end+1} = 'normalize';
        extra_args{end+1} = normalize_values;
    end;
    hold on;
    if iscell(outlinecolor)
        extra_args = [extra_args outlinecolor];
        outlinecolor = color;
    end;
    if isnan(barbottom) || any(orign ~= barbottom) ...
            ... %|| (~isnan(barstart) && any(orign ~= ymin)) % for now, only continue applying this check for deprecated 'barbottom' argument, not new 'barstart' argument
        h(end+1)=plot_bar(origx,orign,outlinecolor,drawn_bin_width,'barbottom',barbottom,'line',extra_args{:});
    else
        h(end+1) = nan;
    end;
end;
