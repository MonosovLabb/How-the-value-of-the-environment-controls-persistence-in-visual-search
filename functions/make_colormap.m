function [cmap] = make_colormap(lo,mid,hi,n)
% [cmap] = make_colormap(lo,mid,hi,n)
%
% code by ESBM, 2011

if nargin < 4
    n = 256;
end;

lo = colorspec_to_rgb(lo);
hi = colorspec_to_rgb(hi);
if isempty(mid)
    mid = lo + hi ./ 2;
else
    mid = colorspec_to_rgb(mid);
end;

nlomid = floor(n/2);
nmidhi = n - nlomid + 1;

rlomid = linspace(lo(1),mid(1),nlomid);
rmidhi = linspace(mid(1),hi(1),nmidhi);

glomid = linspace(lo(2),mid(2),nlomid);
gmidhi = linspace(mid(2),hi(2),nmidhi);

blomid = linspace(lo(3),mid(3),nlomid);
bmidhi = linspace(mid(3),hi(3),nmidhi);

rmidhi = rmidhi(2:end);
gmidhi = gmidhi(2:end);
bmidhi = bmidhi(2:end);

cmap = [[rlomid rmidhi]' [glomid gmidhi]' [blomid bmidhi]'];
