% Y = XE2Y(X,E) - transofrm 3D points X to Y by transofmr E
%
% X ... 3 x n points
% E.r  ... 3 x 3 matrix
% E.t  ... 3 x 1 vector
% E.s  ... 1 x 1 scale
% Y ... 3 x n points:    Y = [E.s*E.r E.t;0 0 0 1]*[X;1]; 

% T. Pajdla, pajdla@cvut.cz, 2017-04-14
function Y = XE2Y(X,E)
if nargin>0
    if size(X,1)<4
        dehom = true;
        X = a2h(X);
    else 
        dehom = false;
    end
    Y = [E.s*E.r E.t;0 0 0 1]*X; 
    if dehom
        Y = h2a(Y);
    end
else % unit tests
    Y = true;
end
