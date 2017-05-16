% f polynomial 1 and 2 formulae, see F2f1f2 for documentation (that is, for
% a link to unexisting documentation.)
y={'(f33)*(f12*f13*f33-f13^2*f32+f22*f23*f33-f23^2*f32)/(f11*f12*f31*f33-f11*f13*f31*f32+f12^2*f32*f33-f12*f13*f32^2+f21*f22*f31*f33-f21*f23*f31*f32+f22^2*f32*f33-f22*f23*f32^2)'
    '(f33)*(f21*f31*f33+f22*f32*f33-f23*f31^2-f23*f32^2)/(f11*f13*f21*f33-f11*f13*f23*f31+f12*f13*f22*f33-f12*f13*f23*f32+f21^2*f23*f33-f21*f23^2*f31+f22^2*f23*f33-f22*f23^2*f32)'
    '(f33)*(f11*f13*f33-f13^2*f31+f21*f23*f33-f23^2*f31)/(f11^2*f31*f33+f11*f12*f32*f33-f11*f13*f31^2-f12*f13*f31*f32+f21^2*f31*f33+f21*f22*f32*f33-f21*f23*f31^2-f22*f23*f31*f32)'
    '(f33)*(f11*f31*f33+f12*f32*f33-f13*f31^2-f13*f32^2)/(f11^2*f13*f33-f11*f13^2*f31+f11*f21*f23*f33+f12^2*f13*f33-f12*f13^2*f32+f12*f22*f23*f33-f13*f21*f23*f31-f13*f22*f23*f32)'};
% f-Ratio formula
y={['f11^2*f13^2*f33-f11*f13^3*f31+2*f11*f13*f21*f23*f33-f11*f13*f23^2*f31+f12^' ...
    '2*f13^2*f33-f12*f13^3*f32+2*f12*f13*f22*f23*f33-f12*f13*f23^2*f32-f13^2*' ...
    'f21*f23*f31-f13^2*f22*f23*f32+f21^2*f23^2*f33-f21*f23^3*f31+f22^2*f23^2*' ...
    'f33-f22*f23^3*f32'],['f11^2*f31^2*f33+2*f11*f12*f31*f32*f33-f11*f13*f31^3-f11*f13*f31*f32^2+f12^' ...
    '2*f32^2*f33-f12*f13*f31^2*f32-f12*f13*f32^3+f21^2*f31^2*f33+2*f21*f22*f31*' ...
    'f32*f33-f21*f23*f31^3-f21*f23*f31*f32^2+f22^2*f32^2*f33-f22*f23*f31^2*f32-' ...
    'f22*f23*f32^3']};
formula=sprintf('%s\n',y{:})
formula=strrep(formula,'f1','F(1,');
formula=strrep(formula,'f2','F(2,');
formula=strrep(formula,'f3','F(3,');
formula=strrep(formula,',1',',1)');
formula=strrep(formula,',2',',2)');
formula=strrep(formula,',3',',3)')