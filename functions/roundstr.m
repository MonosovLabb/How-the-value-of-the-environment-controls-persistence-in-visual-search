function [str] = roundstr(round_digit,num,varargin)
% [str] = roundstr(round_digit,num,...)
%
% convert a number to a string,
%  rounding it to a set digit precision,
%  (e.g. 0.0001)
%  and forcing the string not to be truncated
%  (e.g. force 3.1400 instead of 3.14)
%
%  e.g.:
%  > roundstr(.0001,pi)
%  '3.1416' 
%  > roundstr(.0001,3.14)
%  '3.1400' 
%
% use the num_spaces command to add extra spaces to the start
%  if there are too few characters
% e.g. 
%  > roundstr(1,16,1)
%  '16'
%  > roundstr(1,16,2)
%  '16'
%  > roundstr(1,16,5)
%  '   16'
%  > roundstr(1,16,5,'space','0')
%  '00016'
%  > roundstr(1,16,5,'sign','space','0')
%  '00-16'
%
% code by ESBM, 2012

num_chars = nan;
space_char = ' ';
force_sign = false;

v = 1;
while v <= numel(varargin)
    if ischar(varargin{v})
        switch varargin{v}
            case 'sign'
                force_sign = true;
            case 'space'
                assert(numel(varargin) > v,'"space" argument must be followed by a one-length character argument');
                space_char = varargin{v+1};
                v = v + 1;
            otherwise
                error('unknown argument');
        end;
    else
        assert(isscalar(varargin{v}));
        num_chars = varargin{v};
    end;
    v = v + 1;
end;


num = roundto('digit',round_digit,num);

str = num2str(num);

if round_digit < 1
    n_dec_places = ceil(-log10(round_digit));
    dotpos = find(str=='.',1);
    if isempty(dotpos)
        str = [str '.' repmat('0',1,n_dec_places)];
    else
        n_dec_places_remaining = n_dec_places - (numel(str) - dotpos);
        str = [str repmat('0',1,n_dec_places_remaining)];
    end;
end;
if force_sign && (str(1) ~= '-' && str(1) ~= '+')
    if num < 0 
        str = ['-' str];
    elseif num > 0
        str = ['+' str];
    else
        str = [' ' str];
    end;
end;

if ~isnan(num_chars)
    str = [repmat(space_char,[1 num_chars-numel(str)]) str];
end;
