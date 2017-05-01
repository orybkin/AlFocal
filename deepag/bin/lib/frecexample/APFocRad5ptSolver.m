% C = APFocRad5ptSolver(X,C) - 5pt Minimal Absolute Pose Problem dor principa point at [0,0]
%
% C        = celarray of camera structures, see X2u.m, 'KRCrd' supported 
% X(1:3,:) = 3 x 5 image projections
% X(4:6,:) = 3 x 5 3D points

% T. Pajdla, pajdla@cvut.cz
% 2017-03-06
function C = APFocRad5ptSolver(X,C)

if nargin>0
    if all(size(X)==[6 5])
        C = P5Pfr(X(1:3,:),X(4:6,:));
    end
else % unit tests
    C = true;
end