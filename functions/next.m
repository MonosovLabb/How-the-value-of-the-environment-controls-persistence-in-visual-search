function [v] = next(v,shift,spacer)
% [v] = next(v,shift,spacer)
%
% shifts entries in a vector backwards by 'shift', 
%  filling in empty spaces with the spacer
%  ('spacer' can be a scalar or a vector of length abs(shift))
%
% e.g. next([1 2 3 4],1) == [2 3 4 spacer]
%
% code by ESBM, 2007

if nargin < 2
    shift = 1;
end;
if nargin < 3
    spacer = 0;
end;

v = vshift(v,-shift,spacer);