function [y] = clamp(x,xmin,xmax)
% [y] = clamp(x)
% [y] = clamp(x,xmin,xmax)
%
% clamps the data in x to lie between 'xmin' and 'xmax'
% x, min, and max must be arrays of the same size
%  or scalars
%
% default: xmin = 0, xmax = 1

if nargin == 1
    xmin = 0;
    xmax = 1;
end;

maxdims = max([ndims(x) ndims(xmin) ndims(xmax)]);

sizes = ones(3,maxdims);
sizes(1,1:ndims(x)) = size(x);
sizes(2,1:ndims(xmin)) = size(xmin);
sizes(3,1:ndims(xmax)) = size(xmax);

maxsize = max(sizes);

if isscalar(x) x = x.*ones(maxsize); end;
if isscalar(xmin) xmin = xmin.*ones(maxsize); end;
if isscalar(xmax) xmax = xmax.*ones(maxsize); end;

if ndims(x) ~= maxdims || any(size(x) ~= maxsize) ...
        || ndims(xmin) ~= maxdims || any(size(xmin) ~= maxsize) ...
        || ndims(xmax) ~= maxdims || any(size(xmax) ~= maxsize)
    error('inputs must be scalars or same-sized arrays');
end;

y = max(xmin, min(xmax, x));
