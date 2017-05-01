% E = tabc2tr([x y z a b c]) -  trnslation and rotation angles along XYZ axes to transformation matrix
% 
% [x y z a b c] ... [transpation rotation],  angles a,b,c in degrees
%                    
% E             ... transl([x y z])*rotx(a)*rotx(b)*rotx(c)

% (c) T.Pajdla, www.neovision.cz, Nov 5, 2004
function E = tabc2tr(r)

E = transl(r(1:3))*abc2tr(pi*r(4:6)/180);



