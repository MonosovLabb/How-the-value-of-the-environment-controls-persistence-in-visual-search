function [c] = rowis(x,y)
% function [c] = rowis(x,y)
%  returns a logical column vector,
%  such that c(r) == isequal(x(r,:),y(r,:))

assert(size(x,2)==size(y,2));
if size(x,1)==1 && size(y,1) > 1
    x = repmat(x,size(y,1),1);
end;
if size(x,1) > 1 && size(y,1)==1
    y = repmat(y,size(x,1),1);
end;

c = all(x == y,2);

