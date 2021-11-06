function [v] = vshift(v,shift,spacer)
% [v] = vshift(v,shift,spacer)
%
% shifts entries in a vector forwards (positive shift) or backwards
% (negative shift). 
% if 'spacer' is specified, fills the empty entries with the spacer
%  otherwise, fills them with 'false'
%  ('spacer' can be a scalar or a vector of length abs(shift))
%
% example:
%  > vshift([1 2 3 4],-2)
%  ans = [3 4 0 0]
%  > vshift([1 2 3 4],-2,[99])
%  ans = [3 4 99 99]
%  > vshift([1 2 3 4],-2,[99 43])
%  ans = [3 4 99 43]
%  > vshift([1 2 3 4],+2,[99 43])
%  ans = [99 43 1 2]
%
% code by ESBM, 2008

if nargin < 3 spacer = false; end;

if ~isvector(v) error('v must be a vector'); end;
if ~isscalar(shift) error('shift must be a scalar'); end;
if ~isscalar(spacer) && ~(length(spacer) == abs(shift)) error('spacer must be scalar or a vector with abs(shift) elements'); end;

if isscalar(spacer)
    if isrowvector(v) 
        spacer = spacer*ones(1,abs(shift));
    else
        spacer = spacer*ones(abs(shift),1);
    end;
end;

if shift < 0
    if isrowvector(v) 
        v = v( (1 + abs(shift)) : end);
        v = [v spacer];
    else
        v = v( (1 + abs(shift)) : end);
        v = [v;spacer];
    end;
elseif shift > 0
    if isrowvector(v) 
        v = v( 1 : (end - abs(shift)));
        v = [spacer v];
    else
        v = v( 1 : (end - abs(shift)));
        v = [spacer;v];
    end;
end;

% since the steps above can remove a vector's 'logicalness' 
if all(v == 0 | v == 1)
    v = logical(v);
end;
