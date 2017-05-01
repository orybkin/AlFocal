% y = P4ptP4PfResImPoints(x) - 4pt Minimal Absolute Pose Calibrated up to f - ransac residual evaluation data selection 
%
% x = struct 
%     x.x  = [x1;x2]   ... stacked homogeneous coordinates of image matches
%     x.iK = [iK1;iK2] ... stacked inverses of camera calibration matrices (jssumes triangular iKk)
% y = x.x

% T. Pajdla, pajdla@cvut.cz
% 2017-03-19
function y = P4ptP4PfResImPoints(x)
y = x.x;