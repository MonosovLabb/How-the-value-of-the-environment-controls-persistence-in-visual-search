function [y] = isrowvector(x)
% [y] = isrowvector(x)

y = isvector(x) && size(x,1) == 1;