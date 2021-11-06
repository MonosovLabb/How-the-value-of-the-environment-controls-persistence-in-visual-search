function [h] = setlim(varargin)
% set limits of several axes at once.
%
% [h] = setlim(h,limfxn,lim)
%  for each specified axes handle, call limfxn with lim as its argument.
%  limfxn must be 'xlim', 'ylim', 'xmin','xmax','ymin','ymax', or 'axis'. For example: 
%       setlim(h,'ylim',[0 1]);
%       setlim(h,'axis',[0 .5 999 1001]);
%
% [h] = setlim(h,limfxn,limstyle)
%  set limits to enclose the limits of all plots. 'limstyle' is one of the 
%   following special values:
%    'cur' - use the current limits for each plot
%    'auto' - ask MATLAB to determine the best limits for each plot
%    'tight' - set plot limits to tightly enclose all data points
%     (if there are no data points, displays a warning and does nothing)
%
% [h] = setlim(...,expansion_factor)
%  after doing the above operations, expands the limits by the specified
%  expansion factor.
%   (default = 0).
%
% [h] = setlim(...,expansion_factor,centerpoint,squaremode)
%  after doing the above operations, forces the axes to be symmetric about
%  the specified central point (e.g. 'xlim', centerpoint=0 to force the
%  x-limits to center on 0; or 'axis', centerpoint=[0 5] to force the
%  x-limits to center on 0, and the y-limits to center on 5). If using
%  the 'axis' command, you can specify squaremode as 'square' to
%  additionally force the x and y limits to have the same span
%   (i.e., diff(xlim) == diff(ylim) ). Otherwise, squaremode must be set to
%   'default'.
%   (default: centerpoint = nan, squaremode = 'default').
%
% NOTE: if specify 'h' as 'all', then applies to all axes within the
% current figure.
%
%  examples:
%     setlim(h,@ylim,'tight',.2) 
%       = fit ylimits to data, then expand ymin/ymax by .2*(ymax - ymin)
%
%     setlim(h,@ylim,'tight',[0 .5]) 
%       = fit ylimits to data, then expands only ymax by .5*(ymax - ymin)
%
%     setlim(h,@axis,'tight',[0 .5 .5 0]) 
%       = fit x and y limits to data, then expand xmax by .5*(xmax-xmin),
%         and expand ymin by .5*(ymax-ymin)
%
% code by ESBM, 2011, 2021

% if first argument is not a handle, use the current axis
if numel(varargin) < 1
    error('no input arguments!');
end;
if isempty(varargin{1}) || (~isequal(varargin{1},'all') && ~all(ishandle(varargin{1}(:))))
    varargin = [{gca} varargin];
end;

h = setlim_true(varargin{:});

% -------
function [h] = setlim_true(h,limfxn,lim,expansion_factor,centerpoint,squaremode)

if nargin < 3
    limstyle = 'cur';
end;
if nargin < 4
    expansion_factor = 0;
end;
if nargin < 5
    centerpoint = nan;
end;
if nargin < 6
    squaremode = 'default';
end;


if isempty(h)
    return;
elseif isequal(h,'all')
    h = findobj(gcf,'-property','XLim');
end;

if ischar(lim)
    limstyle = lim;
    lim = [];
    had_only_one_datavalue = [];
    had_no_datavalues = [];
    x_had_only_one_datavalue = nan;
    y_had_only_one_datavalue = nan;
    for hi = 1:numel(h)
        cur_had_no_datavalues = false;
        switch limstyle
            case 'tight'
                h_lines = findobj(h(hi),'-property','XData','-property','YData');
                h_hist = findobj(h(hi),'-property','BinLimits','-property','Values');
                h_text = findobj(h(hi),'units','tight');

                xl = [inf -inf];
                yl = [inf -inf];
                for i = 1:length(h_lines)
                    if ~isequal('off',get(h_lines(i),'xliminclude'))
                        xdat = get(h_lines(i),'XData');
                        xdat = xdat(~isinf(xdat) & ~isnan(xdat)); % ignore infinite and NaN data points
                        if ~isempty(xdat)
                            xl = [min(xl(1),min(xdat(:))) , max(xl(2),max(xdat(:)))];
                        end;
                    end;
                    if ~isequal('off',get(h_lines(i),'yliminclude'))
                        ydat = get(h_lines(i),'YData');
                        ydat = ydat(~isinf(ydat) & ~isnan(ydat)); % ignore infinite and NaN data points
                        if ~isempty(ydat)
                            yl = [min(yl(1),min(ydat(:))) , max(yl(2),max(ydat(:)))];
                        end;
                    end;
                end;
                for i = 1:numel(h_hist)
                    if ~isa(h_hist(i),'matlab.graphics.chart.primitive.Histogram')
                        continue;
                    end
                    if ~isequal('off',get(h_hist(i),'xliminclude'))
                        xdat = get(h_hist(i),'BinLimits');
                        if ~isempty(xdat)
                            xl = [min(xl(1),min(xdat(:))) max(xl(2),max(xdat(:)))];
                        end
                    end
                    if ~isequal('off',get(h_hist(i),'yliminclude'))
                        ydat = get(h_hist(i),'Values');
                        if ~isempty(ydat)
                            yl = [min(yl(1),min(ydat(:))) max(yl(2),max(ydat(:)))];
                        end
                    end
                end
                for i = 1:length(h_text)
                    tpos = get(h_text(i),'Position');
                    tx = tpos(1);
                    ty = tpos(2);
                    if ~isequal('off',get(h_text(i),'xliminclude')) && ~isinf(tx) %ignore infinite data points
                        xl = [min(xl(1),tx) , max(xl(2),tx)];
                    end;
                    if ~isequal('off',get(h_text(i),'yliminclude')) && ~isinf(ty) %ignore infinite data points
                        yl = [min(yl(1),ty) , max(yl(2),ty)];
                    end;
                end;
                % if there was only one data point, remember it
                %  expand the limits by 0.5 on either side of it.
                x_had_only_one_datavalue = nan;
                y_had_only_one_datavalue = nan;
                if xl(2) - xl(1) < eps
                    x_had_only_one_datavalue = xl(1);
                    xl = xl + [-.5 .5];
                end;
                if yl(2) - yl(1) < eps
                    y_had_only_one_datavalue = yl(1);
                    yl = yl + [-.5 .5];
                end;
                

            case 'auto'
                if isequal('auto',get(h(hi),'XLimMode'))
                    xl = xlim(h(hi));
                else
                    curxl = xlim;
                    set(h(hi),'XLimMode','auto');
                    xl = xlim(h(hi));
                    xlim(h(hi),curxl);
                end;

                if isequal('auto',get(h(hi),'YLimMode'))
                    yl = ylim(h(hi));
                else
                    curyl = ylim;
                    set(h(hi),'YLimMode','auto');
                    yl = ylim(h(hi));
                    ylim(h(hi),curyl);
                end;

            case 'cur'
                xl = xlim(h(hi));
                yl = ylim(h(hi));
        end;

        cur_lim = [xl yl];
        cur_had_no_datavalues = any(isinf(xl)) || any(isinf(yl));
        cur_had_only_one_datavalue = [x_had_only_one_datavalue*ones(1,2) y_had_only_one_datavalue*ones(1,2)];
        if hi == 1
            lim = cur_lim;
            had_only_one_datavalue = cur_had_only_one_datavalue;
            had_no_datavalues = cur_had_no_datavalues;
        else
            lim = vertcat(lim,cur_lim);
            had_only_one_datavalue = vertcat(had_only_one_datavalue,cur_had_only_one_datavalue);
            had_no_datavalues = vertcat(had_no_datavalues,cur_had_no_datavalues);
        end;
    end;
    
    
    hadnone = had_no_datavalues==1;
    if all(hadnone)
        warning('setlim couldn''t set tight limits, found no data points!');
        return;
    end;
    lim = lim(~hadnone,:);
    had_only_one_datavalue = had_only_one_datavalue(~hadnone,:);

    % compute overall limits
    newlim = nan*ones(1,size(lim,2));
    for i = 1:size(lim,2)
        curvals = lim(:,i);
                
        % if doing tight fit and there was only one data value,
        % shrink based on that single data value
        % if all plots have the same single data value,
        % shrink to +/- 0.5 around that data value
        hadone = ~isnan(had_only_one_datavalue(:,i));
        curvals(hadone) = had_only_one_datavalue(hadone,i);
        if mod(i,2) == 1
            if all(hadone) && all(curvals - curvals(1) < eps)
                newlim(i) = min(curvals)-0.5;
            else
                newlim(i) = min(curvals);
            end;
        else
            if all(hadone) && all(curvals - curvals(1) < eps)
                newlim(i) = min(curvals)+0.5;
            else
                newlim(i) = max(curvals);
            end;
        end;
    end;
    lim = newlim;
    
    switch limfxn
        case 'xlim'
            lim = lim([1 2]);
        case 'ylim'
            lim = lim([3 4]);
        case 'xmin'
            lim = lim([1]);
        case 'xmax'
            lim = lim([2]);
        case 'ymin'
            lim = lim([3]);
        case 'ymax'
            lim = lim([4]);
        case 'axis'
            lim = lim([1 2 3 4]);
    end;
end;

switch limfxn
    case {'xmin','xmax','ymin','ymax'}
        assert(size(lim,2) == 1);
    case {'xlim','ylim'}
        assert(size(lim,2) == 2);
    case {'axis'}
        assert(size(lim,2) == 4);
end;

% compute how to apply expansion factor to various limits
if numel(expansion_factor) == 1
    expansion_factor = repmat(expansion_factor,1,size(lim,2));
elseif numel(expansion_factor) == 2
    if isequal(limfxn,'axis')
        expansion_factor = repmat(expansion_factor,1,2);
    end;
end;
assert(size(expansion_factor,2) == size(lim,2));

% apply new limits to all axes
for hi = 1:numel(h)
    apply_limfxn(h(hi),limfxn,lim,expansion_factor,centerpoint,squaremode);
end;

% ----
function [curlim] = apply_limfxn(curh,limfxn,limval,expansion_factor,centerpoint,squaremode)

curax = axis(curh);

switch limfxn
    case 'xlim'
        axis(curh,[limval curax([3 4])]);
    case 'ylim'
        axis(curh,[curax([1 2]) limval]);
    case 'xmin'
        axis(curh,[limval curax([2 3 4])]);
    case 'xmax'
        axis(curh,[curax(1) limval curax([3 4])]);
    case 'ymin'
        axis(curh,[curax([1 2]) limval curax(4)]);
    case 'ymax'
        axis(curh,[curax([1 2 3]) limval]);
    case 'axis'
        axis(curh,limval);
end;

curax = axis(curh);
switch limfxn
    case 'xlim'
        axis(curh,[(curax([1 2]) + [-1 1].*expansion_factor.*(curax(2)-curax(1))) curax([3 4])]);
    case 'ylim'
        axis(curh,[curax([1 2]) (curax([3 4]) + [-1 1].*expansion_factor.*(curax(4)-curax(3)))]);
    case 'xmin'
        axis(curh,[(curax(1) + -1.*expansion_factor.*(curax(2)-curax(1))) curax([2 3 4])]);
    case 'xmax'
        axis(curh,[curax(1) (curax(2) + 1.*expansion_factor.*(curax(2)-curax(1))) curax([3 4])]);
    case 'ymin'
        axis(curh,[curax([1 2]) (curax(3) + -1.*expansion_factor.*(curax(4)-curax(3))) curax(4)]);
    case 'ymax'
        axis(curh,[curax([1 2 3]) (curax(4) + 1.*expansion_factor.*(curax(4)-curax(3)))]);
    case 'axis'
        axis(curh,[(curax([1 2]) + [-1 1].*expansion_factor(1:2).*(curax(2)-curax(1))) (curax([3 4]) + [-1 1].*expansion_factor(3:4).*(curax(4)-curax(3)))]);
end;

if ~any(isnan(centerpoint))
    switch limfxn
        case {'xlim','xmin','xmax'}
            assert(numel(centerpoint)==1,'centerpoint must be scalar for x-axis adjustment');
            xl = xlim(curh);
            xlim(curh,centerpoint + [-1 1].*max(abs(xl-centerpoint)));
        case {'ylim','ymin','ymax'}
            assert(numel(centerpoint)==1,'centerpoint must be scalar for y-axis adjustment');
            yl = ylim(curh);
            ylim(curh,centerpoint + [-1 1].*max(abs(yl-centerpoint)));
        case 'axis'
            assert(numel(centerpoint)==1 || numel(centerpoint)==2,'centerpoint must have 1 or 2 entries for axis adjustment');
            if numel(centerpoint) == 1
                centerpoint = [centerpoint centerpoint];
            end;
            xl = xlim(curh);
            yl = ylim(curh);
            max_x_excursion = max(abs(xl-centerpoint(1)));
            max_y_excursion = max(abs(yl-centerpoint(2)));
            switch squaremode
                case 'default'
                    xlim(curh,centerpoint(1) + [-1 1].*max_x_excursion);
                    ylim(curh,centerpoint(2) + [-1 1].*max_y_excursion);
                case 'square'
                    xlim(curh,centerpoint(1) + [-1 1].*max(max_x_excursion,max_y_excursion));
                    ylim(curh,centerpoint(2) + [-1 1].*max(max_x_excursion,max_y_excursion));
                otherwise
                    error('squaremode must be ''default'' or ''square''');
            end;
    
    end;
end;



% ----
function [curlim] = get_lim(curh,limfxn,limval)

curax = axis(curh);

switch limfxn
    case 'xlim'
        curlim = curax([1 2]);
    case 'ylim'
        curlim = curax([3 4]);
    case 'xmin'
        curlim = curax([1]);
    case 'xmax'
        curlim = curax([2]);
    case 'ymin'
        curlim = curax([3]);
    case 'ymax'
        curlim = curax([4]);
    case 'axis'
        curlim = curax;
end;
