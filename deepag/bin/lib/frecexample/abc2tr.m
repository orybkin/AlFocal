% E = abc2tr([a b c]) -  rotation angles along XYZ axes to rotation matrix
% 
% 
% [a b c] ... R = rotx(a)*roty(b)*rotz(c)
% E       ... [R 0; 0 1] rotation matrix

% (c) T.Pajdla, www.neovision.cz, Nov 5, 2004
function R = abc2r(r)

r = r([3 2 1]);
R = rpy2tr(r); 



