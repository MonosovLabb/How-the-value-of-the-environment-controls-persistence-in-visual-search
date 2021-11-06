function [mm] = minmax(x,dim)
% [mm] = minmax(x,dim)
%
%return minimum and maximum of a matrix along dimension dim
%
% if the individual min/max results are row vectors, returns 
%  a them both vertcat'd in a stack:
%   mm = [min_vector ; max_vector]
%
% if the individual min/max results are col vectors, returns 
%  a them both horzcat'd:
%   mm = [min_vector max_vector]
%
% otherwise, adds a dimension onto the matrix, and returns
%  them catt'd in that dimension
%   mm = cat(ndims(min_matrix)+1,min_matrix,max_matrix)
%
% code by ESBM, 2007


if nargin < 2
    min_matrix = min(x);
    max_matrix = max(x);
else
    min_matrix = min(x,[],dim);
    max_matrix = max(x,[],dim);
end;


if isrowvector(min_matrix)
    mm = [min_matrix ; max_matrix];
elseif iscolvector(min_matrix)
    mm = [min_matrix max_matrix];
else
    mm = cat(ndims(min_matrix)+1,min_matrix,max_matrix);
end;
