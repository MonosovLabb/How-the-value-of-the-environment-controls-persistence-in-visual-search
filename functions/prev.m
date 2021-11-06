function [v] = prev(v,shift,spacer)
% [v] = prev(v,shift,spacer)
%
% shifts entries in a vector forwards by 'shift', 
%  filling in empty spaces with the spacer
%  ('spacer' can be a scalar or a vector of length abs(shift))
%
% e.g. prev([1 2 3 4],1) == [spacer 1 2 3]
%
% code by ESBM, 2007

if nargin < 2
    shift = 1;
end;
if nargin < 3
    spacer = 0;
end;

v = vshift(v,shift,spacer);