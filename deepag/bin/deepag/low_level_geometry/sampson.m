function [err,J]=sampson(F,u1,u2)
% calculates the sampson error

if size(u2,1)<3
    u2=[u2; ones(1,size(u2,2))];
    u1=[u1; ones(1,size(u2,2))];
end
F=reshape(F,3,3);
F1=F*u1;
F2=F'*u2;
num=sum(u2.*F1);
den=sqrt((F1(1,:).^2)+(F1(2,:).^2)+(F2(1,:).^2)+(F2(2,:).^2));
err=(num./den);
if nargout>1
    % derivative
    F1(3,:)=0;
    F2(3,:)=0;
    % ((u F v)/den)' = 
    %   (u F v)'*den - num*den'    (u v^T)*den - num*den'
    % = ----------------------- =  ----------------------
    %             den^2                    den^2
    % den' = 1/den * (F1.*v + (u'*F2))
    
    % beware, hacky permutations
    den=permute(den, [1 3 2]);
    num=permute(num, [1 3 2]);
    u2=permute(u2,[1 3 2]);
    u1=permute(u1,[3 1 2]);
    dden = (permute(F1,[1 3 2]).*u1 + u2.*permute(F2, [1 3 2]))./den;
    d=(u2.*u1.*den - num.*dden)./den.^2;
    J=reshape(d,9,[])';
end
end