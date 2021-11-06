function [y] = iscolvector(x)
% [y] = iscolvector(x)

y = isvector(x) && size(x,2) == 1;