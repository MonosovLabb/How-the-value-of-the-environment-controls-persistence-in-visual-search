function y = roundto(method,precis,x)
% y = roundto(method,precis,x)
% example:
%
%%%round to a specific digit:
% roundto('digit',10, 43.0225) == 40
% roundto('digit',1, 43.0225) == 43
% roundto('digit',.1, 43.0225) == 43.0
% roundto('digit',.01, 43.0225) == 43.02
% roundto('digit',.001, 43.0225) == 43.023
%
%%%round to a power of 10:
% roundto('power',1, 43.0225) == 40
% roundto('power',0, 43.0225) == 43
% roundto('power',-1, 43.0225) == 43.0
% roundto('power',-2, 43.0225) == 43.02
% roundto('power',-3, 43.0225) == 43.023
%
%%%round to a power of 10:
% roundto('sigfig',1, 43.0225) == 40
% roundto('sigfig',2, 43.0225) == 43
% roundto('sigfig',3, 43.0225) == 43.0
% roundto('sigfig',4, 43.0225) == 43.02
% roundto('sigfig',5, 43.0225) == 43.023
%
% code by ESBM, 2008

switch method
    case 'digit'
        y = round(x.*(1./precis)).*precis;
    case 'power'
        precis = 10.^precis;
        y = round(x.*(1./precis)).*precis;
    case 'sigfig'
        % get position of first significant figure
        npowers = floor(log10(x));
        precis = 10.^(npowers - precis + 1);
        y = round(x.*(1./precis)).*precis;
    otherwise
        error('unknown method');
end;
