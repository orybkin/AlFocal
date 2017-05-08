function [ F, p, f ] = uu2F_hartley_old(F, u_hartley,v_hartley, f_prior, p_prior, f_min )
% UU2F_HARTLEY solves 7pt problem as described in paper [Reconstruction
% from two views using approximate calibration]
if nargin<6
    f_min=0;
end
% transformation
if size(u_hartley,1)<3
    u_hartley=[u_hartley; ones(1,size(u_hartley,2))];
    v_hartley=[v_hartley; ones(1,size(u_hartley,2))];
end
p_prior=[p_prior{1}; p_prior{2}];
% parameters
niter = 50;
wp=1;
w1=1;
w2=1;
wz1=1;
wz2=1;
wd=1;
%init described in the paper as "Calibrated reconstruction"
F_init=uu2F({u_hartley,v_hartley},{'None' '|F|=0'});
F_init=F;
F=correctF(F_init,f_prior,p_prior);
newF=reshape(F,1,9);
p=p_prior;
f=f_prior;
% options for lsqnonlin
po=optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','off');
fo=optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','off');
Fo=optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','off');
po.Algorithm = 'levenberg-marquardt';
fo.Algorithm = 'levenberg-marquardt';
Fo.Algorithm = 'levenberg-marquardt';

for i=1:niter
    F=newF;    
    p=lsqnonlin(@(p) optimizep(p,p_prior,wp), p,[],[],po);
    f=lsqnonlin(@(f) optimizef(f,f_prior.^2,f_min.^2,w1,w2,wd,wz1,wz2), f.^2,[],[],fo);
    f=sqrt(f);
    newF=lsqnonlin(@(F) sampson(F,u_hartley,v_hartley), F,[],[],Fo);
    newF=reshape(correctF(reshape(newF,3,3),f,p),1,9);
    if norm(newF-F)<10^-6
        break
    end
    %iter(i)=
end
F=newF;
p={p(1:2) p(3:4)};
end

function [err, J]=optimizep(p,p_prior,wp)
err=wp*(p-p_prior);
if nargout>1
    J=diag(repmat(wp,size(p,1),1));
end
end

function [err, J]=optimizef(f,f_prior,f_min,w1,w2,wd,wz1,wz2)
err=[w1*(f(1)-f_prior(1))
      w2*(f(2)-f_prior(2))    
      wd*(f(1)-f(2))
      wz1*(f(1)-f_min)*(f(1)<f_min)
      wz2*(f(2)-f_min)*(f(2)<f_min)];
if nargout>1
    J=[w1 0
       0 w2 
       wd -wd
       wz1*(f(1)<f_min) 0
       0 wz2*(f(2)<f_min)];
end
end



