function [rgb] = colorspec_to_rgb(cs)
% [rgb] = colorspec_to_rgb(cs)
%  convert a color specification to an RGB triplet
%  if doesn't match a full color name like 'magenta',
%  strips down to only the characters 'ymcrgbwk' and tries
%  again. So, you can pass an argument to the plot function,
%  like ':go', and it'll correctly extract the green triplet [0 1 0]
%
% Code by ESBM, 2007

if ischar(cs)
    cs = lower(cs);
    switch cs
        case {'y','yellow'}
            rgb = [1 1 0];
        case {'m','magenta'}
            rgb = [1 0 1];
        case {'c','cyan'}
            rgb = [0 1 1];
        case {'r','red'}
            rgb = [1 0 0];
        case {'g','green'}
            rgb = [0 1 0];
        case {'b','blue'}
            rgb = [0 0 1];
        case {'w','white'}
            rgb = [1 1 1];
        case {'k','black'}
            rgb = [0 0 0];
        otherwise
            % if there are extra characters, strip them and try again
            strippedcs = cs(regexp(cs,'[ymcrgbwk]'));
            if strcmp(strippedcs,cs)
                error('invalid colorspec!');
            else
                rgb = colorspec_to_rgb(strippedcs);
            end;
    end;
elseif ndims(cs) == 2 && all(size(cs) == [1 3])
    rgb = cs;
else
    error('invalid colorspec!');
end;