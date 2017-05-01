% e = xy2HSolverRes(H,X) - 2D Homography Solver Residuals
%
% H        = homography matrix celarray
% X(1:2,:) = 2 x n points in plane 1
% X(3:4,:) = 2 x n points in plane 2
% e        = residuals (xy2HRes for 'y-y')

% T. Pajdla, 2016-09-01
function e = xy2HSolverRes(H,X,~,~)
if nargin>0
    % e = sum(sum((X(3:4,:)-h2a(H{1}*a2h(X(1:2,:)))).^2));
    e = xy2Hres(H{1},a2h(X(1:2,:)),a2h(X(3:4,:)),'y-y')'; % must be a row vector for each model
else % unit tests
    e = true;
end
