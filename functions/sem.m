function [x_sem] = sem(x,dim)
% [x_sem] = sem(x,dim)
%calculate SEM of data
% along dimension 'dim'

if nargin < 2
    x_sem = std(x,0) ./ sqrt(sum(ones(size(x))));
else
    x_sem = std(x,0,dim) ./ sqrt(sum(ones(size(x)),dim));
end;

