function [icolor] = interpcolor(color1,color2,alpha)
% [icolor] = interpcolor(color1,color2,alpha)
% interpolate between two colors, simply using:
%  icolor = color1.*((1-alpha).*[1 1 1]) + color2.*((alpha).*[1 1 1])
% default: alpha = .5
% 
% code by ESBM, 2011

if nargin < 2
    color2 = color1;
end;
if nargin < 3
    alpha = .5;
end;


icolor = colorspec_to_rgb(color1).*((1-alpha).*[1 1 1]) + colorspec_to_rgb(color2).*((alpha).*[1 1 1]);
