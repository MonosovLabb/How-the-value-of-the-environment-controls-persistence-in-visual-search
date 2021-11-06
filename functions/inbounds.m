function [y] = inbounds(x,bounds,edgetype)
% [y] = inbounds(x,bounds)
%
% x = variable to be tested
% bounds = [lower upper] bounds. Can be specified as a 1x2 matrix, or if
%          x is a Nx1 column vector, can be a N x 2 matrix, where each row
%          of bounds is applied to the corresponding entry in x.
% edgetype = 2-length string composed of:
%  '[' ']' for left or right closed bounds (i.e. <=)
%  '(' ')' for left or right open bounds (i.e. <)
% (default: '[]')
% 
% y = same size as x. Each elemnt is true if the corresponding
%  element of x is in the bounds 'bounds'
% clamps the data in x to lie between 'xmin' and 'xmax'
%
% example:
%  > inbounds([1 2 3 4 5 6],[2 5],'[)')
%  ans = [0 1 1 1 0 0]
% 
% code by ESBM, 2018

if nargin < 3
    edgetype = '[]';
end

if iscolumn(x)
    assert(isequal(size(bounds),[1 2]) || isequal(size(bounds),[size(x,1) 2]),'if x is a column vector, bounds must be 1x2 or Nx2');
else
    assert(isequal(size(bounds),[1 2]),'if x is not a column vector, bounds must be 1x2');
end
assert(all(bounds(:,1) <= bounds(:,end)) && ischar(edgetype) && numel(edgetype)==2,'bounds(:,1) must be <= bounds(:,end) and edgetype must be a 2-length character array');

switch edgetype
    case '[]'
        y = bounds(:,1) <= x & x <= bounds(:,2);
    case '()'
        y = bounds(:,1) <  x & x <  bounds(:,2);
    case '(]'
        y = bounds(:,1) <  x & x <= bounds(:,2);
    case '[)'
        y = bounds(:,1) <= x & x <  bounds(:,2);
    otherwise
        error('edgetype must be [] or () or (] or [)');
end