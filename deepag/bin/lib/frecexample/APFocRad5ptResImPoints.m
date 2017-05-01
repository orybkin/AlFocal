% y = APFocRad5ptResImPoints(x) - 5pt Minimal Absolute Pose Problem dor principa point at [0,0]
%
% x = struct 
%     x.u = 3 x n homogeneous cordinates of imahe points
%     x.X = 3 x n cordinates of eD points
%     x.C = camera paramteres
% y = [x.u; x.X]

% T. Pajdla, pajdla@cvut.cz
% 2017-03-19
function y = APFocRad5ptResImPoints(x)
K = [1 0 x.C.K(1,3); 0 1 x.C.K(2,3); 0 0 1]; 
y = [K\a2h(x.u); x.X];
