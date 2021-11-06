function [J,clim] = colormapify(I,clim,lo,mid,hi,nancolor)
% [J] = colormapify(I,clim,cmap)
% given an (r x c) image, an intensity range to
%  spread colors over, and a colormap,
%  create the resulting (r x c x 3) rgb image
%
% if a pixel is higher than clim(2), assigns
%  it the 'hi' color; if a pixel is lower than clim(1),
%  assigns it the 'lo' color.
%
% [J] = colormapify(I,[],cmap)
%  defaults to use colormap spanning the full range of data
%
% [J] = colormapify(I,clim,lo,mid,hi)
% creates the colormap using the function make_colormap
%
% [J] = colormapify(I,clim,lo,mid,hi,outrange)
% creates the colormap using the function make_colormap
%  and sets NaN pixels to the color 'nancolor'
%  (default: the 'middle' color of the colormap)
%
% Code by ESBM, 2008

% size of the colormap
n = 256;

if isempty(clim)
    clim = [min(I(:)) max(I(:))];
end;

if nargin < 4
    cmap = lo;
else
    cmap = make_colormap(lo,mid,hi);
end;
if nargin < 6
    nancolor = cmap(round(size(cmap,1)./2),:);
else
    nancolor = colorspec_to_rgb(nancolor);
end;

cbins = linspace(clim(1),clim(2),n+1);
[jnk,ids] = histc(I(:),cbins);

ids(ids==n+1) = n;
ids(I(:) < clim(1)) = 1;
ids(I(:) > clim(2)) = n;
ids = reshape(ids,size(I));

J = zeros(size(I,1),size(I,2),3);
for k = 1:3
    for r = 1:size(I,1)
        goodids = 1 <= ids(r,:) & ids(r,:) <= n;
        J(r,goodids,k) = cmap(ids(r,goodids),k);
        J(r,~goodids,k) = nancolor(1,k);
        % J(r,:,k) = cmap(ids(r,:),k);
    end;
end;
