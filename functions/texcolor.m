function [str] = texcolor(rgb)
% [str] = texcolor(rgb)
% produce a TeX string used to change
% the text color. e.g.
% coloredtext = [texcolor([0 .5 1]) '{' oldtext '}'];
%
% code by ESBM, 2007


rgb = colorspec_to_rgb(rgb);
str = ['\color[rgb]{'  num2str(rgb(1)) ' ' num2str(rgb(2)) ' ' num2str(rgb(3))  '}'];
