% H = xy2HSolver(X[,H]) - 2D Homography Solver
%
% H        = homography matrix celarray
% X(1:2,:) = 2 x n points in plane 1
% X(3:4,:) = 2 x n points in plane 2

% T. Pajdla, 2016-09-01
function H = xy2HSolver(X,~)
if nargin>0
    H = {xy2H(a2h(X(1:2,:)),a2h(X(3:4,:)),{'DLT','y-y'})};
else % unit tests
    H = true;
end
