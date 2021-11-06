function [n] = nans(varargin)
% [n] = nans(...)
% returns a matrix of NaNs
% nans(...) is the same as nan*ones(...)

n = nan*ones(varargin{:});
