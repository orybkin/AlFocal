% y = P4ptP4PfFitImPoints(x) - 4pt Minimal Absolute Camera Pose - ransac calibrated up to f fit data computation
%
% x = struct 
%     x.x  = homogeneous coordinates of image matches
%     x.iK = inverse of camera calibration matrix
% y = iK*x.x

% T. Pajdla, pajdla@cvut.cz
% 2017-03-19
function y = P4ptP4PfFitImPoints(x)
y = [x.iK*x.x(1:3,:)/x.iK(3,3); x.x(4:end,:)];