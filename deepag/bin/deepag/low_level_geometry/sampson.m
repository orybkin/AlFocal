function [err,J]=sampson(F,u,v)
% calculates the sampson error

if size(u,1)<3
    u=[u; ones(1,size(u,2))];
    v=[v; ones(1,size(u,2))];
end
F=reshape(F,3,3);
F1=F*v;
F2=F'*u;
num=sum(u.*F1);
den=sqrt((F1(1,:).^2)+(F1(2,:).^2)+(F2(1,:).^2)+(F2(2,:).^2));
err=(num./den);
if nargout>1
    F1(3,:)=0;
    F2(3,:)=0;
    % ((u F v)/den' = 
    %   (u F v)'*den - num*den'    (u v^T)*den - num*den'
    % = ----------------------- =  ----------------------
    %             den^2                    den^2
    % den' = 1/den * (F1.*v + (u'*F2))
    
    % beware, hacky permutations
    den=permute(den, [1 3 2]);
    num=permute(num, [1 3 2]);
    u=permute(u,[1 3 2]);
    v=permute(v,[3 1 2]);
    dden = (permute(F1,[1 3 2]).*v + u.*permute(F2, [1 3 2]))./den;
    d=(u.*v.*den - num.*dden)./den.^2;
    J=reshape(d,9,[])';
end
end