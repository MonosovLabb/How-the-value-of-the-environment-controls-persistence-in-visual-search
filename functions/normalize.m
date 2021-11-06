function [x] = normalize(x,n)
% [x] = normalize(x,n)
%
%normalize x so that its n-norm is equal to 1.
% x must be a vector
%
%Example:
% if n = 1 and x is non-negative, 
% then x becomes a probability distribution, 
% i.e. with sum(x) = 1
%Example:
% if n = 2,
% then x becomes a unit vector,
% i.e. with euclidean distance = 1
%
% (default: n = 2)
%
% code by ESBM, 2007

if ~isvector(x)
    error('x must be a vector');
end;

x = double(x);

if nargin < 2 
    x = x ./ norm(x);
else
    x = x ./ norm(x,n);
end;
