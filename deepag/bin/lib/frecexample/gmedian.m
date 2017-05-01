% M = gmedian(X) - geometric median of kD points X
% 
% X ... k x n point coordinates
% M ... k x 1 geometric median of X
%             enumerative implementation, see the following for other approaches
%             https://en.wikipedia.org/wiki/Geometric_median#CITEREFCohenLeeMillerPachocki2016
%             https://arxiv.org/pdf/math/0702005.pdf 
%             https://arxiv.org/pdf/1606.05225.pdf

% (c) T.Pajdla, pajdla@cvut.cz, 2017-04-15
function M = gmedian(X)
if nargin>0
    d = XdX(X); % pairwise distances
    d = sum(d); % sum of distances
    i = argmin(d); % the optimal point index
    M = X(:,i); % the optimal point
else % unit tests 
    X = [0 1 1 0 0.5
         0 0 1 1 0.5
         1 1 1 1 1]; 
    M = gmedian(X);
    M = ~any(M-[0.5;0.5;1]);
end