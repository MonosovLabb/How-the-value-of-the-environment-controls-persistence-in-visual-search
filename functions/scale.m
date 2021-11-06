function [h] = scale(varargin)
% [h] = scale(varargin)
% for rescaling figures/axes to meet
%  specified requirements
%
% arguments:
%  handle: 
%   [matrix of figure handles]
%   'fig',[matrix of figure handles]
%    (in which case scales each figure independently)
%
%   [matrix of axes handles from the same figure]
%   'ax',[matrix of axes handles from the same figure]
%    (in which case, scales the bounding box of the collection of axes)
%
%   'ax',[matrix of figure handles]
%    (selects all the axes in each figure, and scales their bounding box)
%
%   'eachax',[matrix of figure handles]
%    (scales each axes independently of the others)
%
%  scale reference point: any one of nine points
%   horizontal position: 'l'eft, 'c'enter, 'r'ight
%   vertical position: 'b'ottom, 'm'iddle, 't'op
%   (can also specify as -1/0/1 instead of 'l'/'c'/'r' or 'b'/'m'/'t')
%  
%   (or if in 'fig' mode, can be 'center', in which case centers the
%    figures on the screen and then returns)
%   (NOTE: this argument is optional. If not specified, defaults to 'ct')
%
%  aspect mode:
%   'only': make scaling changes to horz/vert dimensions independently
%   'lock': make same scaling changes to horz/vert dimensions (preserve aspect ratio)
%    (overridden if explicitly give diff scaling factors in each dimension)
%   (NOTE: this argument is optional. If not specified, defaults to 'only')
%
%  scale mode: (when scaling axes only)
%   'pos'ition: scale the 'position' property of the axes (the default)
%   'out'erposition: scale the 'outer position' property of the axes
%   'tight'inset: scale the 'tightinset' property. Cannot do this
%     with a direct command. Instead, does it by scaling 'position'
%     incrementally.
%
%  contraint: (multiple arguments)  
%   'hlen',[new length],[units]
%   'vlen',[new length],[units]
%   'len' (or 'hvlen'),[new length],[units] (makes it square)
%   'len' (or 'hvlen'),[new h length, new v length],[units]
%     (if units is unspecified, uses current figure/axis units)
%   'hscale',[scale factor]
%   'vscale',[scale factor]
%   'scale' (or 'hvscale'),[scale factor]
%   'scale' (or 'hvscale'),[h scale factor, v scale factor]
%
%    units can be: in(ches), cm, mm, cen(timeters), norm(alized), px (pixels), pt, point
%
% examples:
%  scale figure from center, keeping aspect ratio, to have horz. length = 129 mm
%   scale(fig_handle,'cm','lock','hlen',129,'mm')
%  scale figure from lower side, keeping aspect ratio, to have vertical. length = 129 mm
%   scale(fig_handle,'cm','lock','hlen',129,'mm')
%
%
% NOTE: for axes, may have to reprogram so that 
%  if set 'Position' will also set ActivePositionProperty = 'Position',
%   (same for TightInset?)
%  and if set 'OuterPosition' will also set it to 'OuterPosition'.
%
% code by ESBM, 2012, 2014, 2021

h = [];

is_fig_handle = @(h) get(h,'parent')==0;
is_ax_handle = @(h) ~isempty(findobj(h,'-depth',0,'-property','xlim'));

mode = 'fig';
fig_handle = [];
ax_handles = [];


scalerefpoint = [0 0];
scaleaspect = 'lock';
scalemode = 'pos';
scaleconstraint = {};

% --- read handle mode ---
v = 1;
if ischar(varargin{v})
    switch varargin{v}
        case {'fig','figure'}
            mode = 'fig';
            v = v + 1;
        case {'ax','axes'}
            mode = 'ax';
            v = v + 1;
        case {'eachax','axes'}
            mode = 'eachax';
            v = v + 1;
        otherwise
            error('unknown handle type');
    end;
else
    if ~isempty(varargin{v}) && is_fig_handle(varargin{v}(1))
        mode = 'fig';
    elseif ~isempty(varargin{v}) && is_ax_handle(varargin{v}(1))
        mode = 'ax';
    else
        error('unknown handle type');
    end;
end;

% --- read handles ---
recursive_args = varargin(v+1:end);
switch mode
    case 'fig'
        for i = 1:numel(varargin{v})
            assert(is_fig_handle(varargin{v}(i)),'expected figure handles');
        end;
        fig_handle = varargin{v};
        
        if isempty(fig_handle)
            % no handles? then nothing to do
            return;
        elseif numel(fig_handle) > 1
            % call recursively on each figure handle
            for fh = 1:numel(fig_handle)
                h = [h scale(fig_handle(fh),recursive_args{:})];
            end;
            return;
        end;
        
        % otherwise, we have just a single figure
        %  do not call recursively. We'll handle it in this function!
        
    case {'eachax','ax'}
        % ignore NaN handles
        h = h(~isnan(h));
        
        isvalid = false(size(varargin{v}));
        isfh = false(size(varargin{v}));
        isah = false(size(varargin{v}));
        for i = 1:numel(varargin{v});
            isvalid(i) = isa(varargin{v}(i),'matlab.graphics.axis.Axes') || isa(varargin{v}(i),'matlab.ui.Figure') || ~isnan(varargin{v}(i));
            if isvalid(i)
                isfh(i) = is_fig_handle(varargin{v}(i));
                isah(i) = is_ax_handle(varargin{v}(i));
            end;
        end;
        
        if all(~isvalid)
            disp(['scale received all NaN handles, doing nothing']);
            return;
        elseif all(isfh(isvalid))
            % call recursively on each figure handle
            ax_fig_handles = varargin{v};
            for afh = 1:numel(ax_fig_handles)
                % find child axes of the figure
                children = get(ax_fig_handles(afh),'Children');
                childax = false(size(children));
                for i = 1:numel(children)
                    childax(i) = is_ax_handle(children(i));
                end;
                children = children(childax);
                % call recursively on the child axes
                h = [h scale(mode,children,recursive_args{:})];
            end;
            % done!
            return;
            
        elseif all(isah(isvalid))
            % if axes handles are from more than one figure,
            %  call recursively on axis handles from the same figure
            ax_handles = varargin{v};
            ax_fig = nans(size(ax_handles));
            for ah = 1:numel(ax_handles)
                ax_fig(ah) = get(ax_handles(ah),'Parent');
                assert(is_fig_handle(ax_fig(ah)),'expected parent of axes to be a figure');
            end;
            unique_figs = unique(ax_fig);
            if numel(unique_figs) > 1
                for uf = 1:numel(unique_figs)
                    h = [h scale(mode,ax_handles(ax_fig==unique_figs(uf)),recursive_args{:})];
                end;
                return;
            elseif isequal(mode,'eachax')
                % if in 'eachax' mode,
                %  call scale recursively on each axis independently
                for i = 1:numel(ax_handles)
                    h = [h scale('ax',ax_handles(i),recursive_args{:})];
                end;
                return;
            end;
            
            % otherwise, axes handles are all from the same figure
            %  do not call recursively. We'll handle them in this function!
            
        else
            error('expected all figure handles or all axes handles');
        end;
        
        % no handles? then nothing to do
        if isempty(ax_handles)
            return;
        end;
end;

assert(isequal(mode,'fig') && ~isempty(fig_handle) && isempty(ax_handles) ...
    || isequal(mode,'ax') && isempty(fig_handle) && ~isempty(ax_handles), 'expected only fig_handle in ''fig'' mode or only ax_handles in ''ax'' mode');
if isequal(mode,'ax')
    fig_handle = get(ax_handles(1),'parent');
    assert(is_fig_handle(fig_handle),'expected parent of an axis to be a figure');
end;

v = v + 1;

% --- read scale reference point ---
if ~isempty(varargin{v})
    if ischar(varargin{v})
        
        % center figure on screen
        if isequal(mode,'fig') && isequal(varargin{v},'center')
            center_figure_on_screen(fig_handle);
            h = fig_handle;
            return;
        end;
        

        if numel(varargin{v})==1
            switch lower(varargin{v}(1))
                case 'l'
                    scalerefpoint = [-1 0];
                case 'c'
                    scalerefpoint = [0 0];
                case 'r'
                    scalerefpoint = [1 0];
                case 'b'
                    scalerefpoint = [0 -1];
                case 'm'
                    scalerefpoint = [0 0];
                case 't'
                    scalerefpoint = [0 1];
                otherwise
                    error('unexpected input to scalerefpoint');
            end;
        elseif numel(varargin{v})==2
            switch lower(varargin{v}(1))
                case 'l'
                    scalerefpoint(1) = -1;
                case 'c'
                    scalerefpoint(1) = 0;
                case 'r'
                    scalerefpoint(1) = 1;
                otherwise
                    error('unknown first input to scalerefpoint');
            end;
            switch lower(varargin{v}(2))
                case 'b'
                    scalerefpoint(2) = -1;
                case 'm'
                    scalerefpoint(2) = 0;
                case 't'
                    scalerefpoint(2) = 1;
                otherwise
                    error('unknown second input to scalerefpoint');
            end;
        else
            % default! scalerefpoint = 'cm'
            scalerefpoint = [0 1];
            v = v - 1;
        end;
    else
        % default! scalerefpoint = 'cm'
        scalerefpoint = [0 1];
        v = v - 1;
    end;
end;

v = v + 1;

% --- read aspect mode ---
if ~isempty(varargin{v})
    if ischar(varargin{v})
        scaleaspect = lower(varargin{v});
        switch scaleaspect
            case {'only','lock'}
            otherwise
                % default! scaleaspect = 'only'
                scaleaspect = 'only';
                v = v - 1;
        end;
    else
        % default! scaleaspect = 'only'
        scaleaspect = 'only';
        v = v - 1;
    end;
end;

v = v + 1;

% --- read scalemode (in 'axes' mode) ---
if isequal(mode,'ax')
    if ~isempty(varargin{v})
        scalemode = lower(varargin{v});
        switch scalemode
            case {'pos','out'}
                v = v + 1;
            otherwise
                scalemode = 'pos';
        end;
    end;
end;


% --- read constraint ---
scaleconstraint = varargin(v:end);

switch mode
    % --- figure scaling ---
    case 'fig'
        h = fig_handle;
        
        if numel(scaleconstraint)<2
            return;
        end;
        
        sctype = lower(scaleconstraint{1});
        val = scaleconstraint{2};
        assert(~isempty(sctype),'need a non-empty scale constraint type');
        assert(~isempty(val),'need a non-empty scale constraint value');
        
        % temporarily change figure units
        old_units = get(fig_handle,'units');
        if numel(scaleconstraint)>=3
            new_units = scaleconstraint{3};
            switch new_units
                case 'cm'
                    new_units = 'centimeters';
                case 'in'
                    new_units = 'inches';
                case 'pt'
                    new_units = 'point';
                case 'px'
                    new_units = 'pixels';
                case 'mm'
                    new_units = 'centimeters';
                    switch sctype
                        case {'hlen','vlen','hvlen'}
                            val = val ./ 10;
                    end;
            end;
        else
            new_units = old_units;
        end;
        set(fig_handle,'units',new_units);
        
        
        old_lbwh = get(fig_handle,'pos');
        old_lrbt = lbwh_to_lrbt(old_lbwh);
        [refpoint,scalefactor] = convert_constraints_to_scaling(old_lrbt,scalerefpoint,scaleaspect,sctype,val);
        new_lrbt = apply_scaling(old_lrbt,refpoint,scalefactor);
        
        % rescale figure
        new_lbwh = lrbt_to_lbwh(new_lrbt);
        
        %2014-02-01 seems like matlab has problems scaling the figure when
        % the x-coordinate is a negative number. Even though this is
        %  still possibly a valid position (e.g. on a diff display of
        %  a multi-display monitor). Perhaps the true problem is that it
        %  can't scale a figure if the figure is stretched across
        %  two displays?
        if new_lbwh(1) < 1
            new_lbwh(1) = 1;
        end;
        
        set(fig_handle,'pos',new_lbwh);
        waittimer = tic;
        while toc(waittimer) < 0.1
            % wait 100 ms for the 'position change' to take effect
            %  before checking it. Oddly, this seems to be necessary.
            % If you resize the figure to be very large, its height/width
            %  get clamped to a maximum value (how is it decided??)
            continue;
        end;
        real_new_lbwh = get(fig_handle,'pos');
        if ~isequal(new_lbwh,real_new_lbwh)
            warning('figure position not changed correctly - centering figure on screen. Figure position may become unstable.');
            center_figure_on_screen(fig_handle);
        end;
        
        % reset to old units
        set(fig_handle,'units',old_units);
        
        
    case 'ax'
        h = ax_handles;
        
        if numel(scaleconstraint)<2
            return;
        end;
        
        sctype = lower(scaleconstraint{1});
        val = scaleconstraint{2};
        assert(~isempty(sctype),'need a non-empty scale constraint type');
        assert(~isempty(val),'need a non-empty scale constraint value');
        
        % temporarily change axis units, if necessary
        old_units = get(ax_handles,'units');
        if ~iscell(old_units)
            old_units = {old_units};
        end;
        if numel(scaleconstraint)>=3
            switch scaleconstraint{3}
                case 'cm'
                    new_units = 'centimeters';
                case 'in'
                    new_units = 'inches';
                case 'pt'
                    new_units = 'point';
                case 'px'
                    new_units = 'pixels';
                case 'mm'
                    new_units = 'centimeters';
                    switch sctype
                        case {'hlen','vlen','hvlen'}
                            val = val ./ 10;
                    end;
            end;
            new_units = repmat({new_units},size(old_units));
        else
            new_units = old_units;
        end;
        for i = 1:numel(ax_handles)
            set(ax_handles(i),'units',new_units{i});
        end;
        
        % figure out desired new bounding box (for all axes combined)
        old_lrbt = get_ax_bounding_box(ax_handles,scalemode);
        [refpoint,scalefactor] = convert_constraints_to_scaling(old_lrbt,scalerefpoint,scaleaspect,sctype,val);
        new_lrbt = apply_scaling(old_lrbt,refpoint,scalefactor);
        
        % get desired new bounding box for each individual axes,
        %  and then change its position to fit that bounding box
        ax_old_lbwh = nans(numel(ax_handles),4);
        ax_old_lrbt = nans(numel(ax_handles),4);
        ax_new_lrbt = nans(numel(ax_handles),4);
        ax_new_lbwh = nans(numel(ax_handles),4);
        for i = 1:numel(ax_handles)
            ax_old_lbwh(i,:) = get(ax_handles(i),scalemode);
            ax_old_lrbt(i,:) = lbwh_to_lrbt(ax_old_lbwh(i,:));
            ax_new_lrbt(i,:) = apply_scaling(ax_old_lrbt(i,:),refpoint,scalefactor);
            ax_new_lbwh(i,:) = lrbt_to_lbwh(ax_new_lrbt(i,:));
        end;
        
        switch scalemode
            case {'pos','out'}
                % rescale the 'position' or 'outerposition' properties
                %  directly
                for i = 1:numel(ax_handles)
                    set(ax_handles(i),scalemode,ax_new_lbwh(i,:));
                end;
            
            otherwise
                error('unknown scalemode');
        end;
        
        % return axes to old units
        for i = 1:numel(ax_handles)
            set(ax_handles(i),'units',old_units{i});
        end;
    otherwise
        error('unknown mode, expected ''fig'' or ''ax''');
end;


% ----
function [lrbt] = lbwh_to_lrbt(lbwh)
lrbt = [lbwh(1) lbwh(1)+lbwh(3) lbwh(2) lbwh(2)+lbwh(4)];

% ----
function [lbwh] = lrbt_to_lbwh(lrbt)
lbwh = [lrbt(1) lrbt(3) diff(lrbt(1:2)) diff(lrbt(3:4))];


% ----
function [] = center_figure_on_screen(fig_handle)
old_units = get(fig_handle,'units');

monitor_pos = get(0,'MonitorPosition');
set(fig_handle,'units','pixels');
curpos = get(fig_handle,'pos');
set(fig_handle,'pos',[((monitor_pos(3)./2)-(curpos(3)./2)) ((monitor_pos(4)./2)-(curpos(4)./2)) curpos(3:4)]);

set(fig_handle,'units',old_units);

% ----
function [lrbt] = get_ax_bounding_box(ax_handles,scalemode)

% left right bottom top
lrbt = [nan nan nan nan];

% left bottom width height
lbwh_ax = nans(numel(ax_handles),4);
for ah = 1:numel(ax_handles)
    lbwh_ax(ah,:) = get(ax_handles(ah),scalemode); %scalemode is 'pos'ition | 'out'erposition | 'tight'inset
end;

% convert from lbwh to lrbt
lrbt_ax = [lbwh_ax(:,1) lbwh_ax(:,1)+lbwh_ax(:,3) lbwh_ax(:,2) lbwh_ax(:,2)+lbwh_ax(:,4)];

% find bounding box
lrbt = [min(lrbt_ax(:,1)) max(lrbt_ax(:,2)) min(lrbt_ax(:,3)) max(lrbt_ax(:,4))];


% ----
function [refpoint,scalefactor] = convert_constraints_to_scaling(old_lrbt,scalerefpoint,scaleaspect,sctype,val)
% convert various constraints to a standard transformation:
%  scaling with reference to a specified reference point in figure coords.

% find true reference point
refpoint = [nan nan];
switch scalerefpoint(1)
    case -1
        refpoint(1) = old_lrbt(1);
    case 0
        refpoint(1) = mean(old_lrbt(1:2));
    case 1
        refpoint(1) = old_lrbt(2);
    otherwise
        error('unknown scalerefpoint(1)');
end;
switch scalerefpoint(2)
    case -1
        refpoint(2) = old_lrbt(3);
    case 0
        refpoint(2) = mean(old_lrbt(3:4));
    case 1
        refpoint(2) = old_lrbt(4);
    otherwise
        error('unknown scalerefpoint(2)');
end;

% convert length changes to scalings
switch sctype
    case 'hlen'
        assert(isscalar(val));
        val = (val ./ diff(old_lrbt(1:2)));
        sctype = 'hscale';
    case 'vlen'
        assert(isscalar(val));
        val = (val ./ diff(old_lrbt(3:4)));
        sctype = 'vscale';
    case {'len','hvlen'}
        if isscalar(val)
            val = [val val];
        end;
        assert(numel(val)==2);
        val(1) = (val(1) ./ diff(old_lrbt(1:2)));
        val(2) = (val(2) ./ diff(old_lrbt(3:4)));
        sctype = 'hvscale';
end;

% convert horz/vert scalings to combined horz/vert scaling
switch sctype
    case 'hscale'
        assert(isscalar(val));
        if isequal(scaleaspect,'lock')
            val = [val val];
        else
            val = [val 1];
        end;
        sctype = 'hvscale';
    case 'vscale'
        assert(isscalar(val));
        if isequal(scaleaspect,'lock')
            val = [val val];
        else
            val = [1 val];
        end;
        sctype = 'hvscale';
    case 'scale'
        sctype = 'hvscale';
end;

assert(isequal(sctype,'hvscale'));
if isscalar(val)
    val = [val val];
end;
assert(numel(val)==2);

% x,y scale factors
scalefactor = [val(1) val(2)];


% ----
function [new_lrbt] = apply_scaling(old_lrbt,refpoint,scalefactor)
% scale a bounding box with respect to a fixed (x,y) reference point
%  and a specified scale factor

new_lrbt = old_lrbt;
new_lrbt(1) = refpoint(1) + (old_lrbt(1) - refpoint(1))*scalefactor(1);
new_lrbt(2) = refpoint(1) + (old_lrbt(2) - refpoint(1))*scalefactor(1);
new_lrbt(3) = refpoint(2) + (old_lrbt(3) - refpoint(2))*scalefactor(2);
new_lrbt(4) = refpoint(2) + (old_lrbt(4) - refpoint(2))*scalefactor(2);

