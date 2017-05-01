% function s = pp(x) - projectors s{i} = x(:,i)*x(:,i)'
% 
% x ...  3 x n matrix [a b c; ...]'
% s ... { x(:,i)*x(:,i)' ... }

% T. Pajdla, pajdla@cvut.cz
% (c) 2017-03-05
function s = pp(x)
if size(x,2)==1
    s{1} = x*x';   
else
    for i=1:size(x,2)
        s{i} = x(:,i)*x(:,i)';   
    end
end