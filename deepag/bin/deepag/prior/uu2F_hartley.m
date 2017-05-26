function [ F, p, f, output ] = uu2F_hartley(F, u,v, f_prior, p_prior, f_min, f_ratio )
% UU2F_HARTLEY solves 7pt problem as described in paper [Reconstruction
% from two views using approximate calibration]
if nargin<6
    f_min=0;
end
if nargin<7
    f_ratio=1;
end
% transformation
if size(u,1)<3
    u=[u; ones(1,size(u,2))];
    v=[v; ones(1,size(u,2))];
end
p_prior=[p_prior{1}; p_prior{2}];
% parameters
wF=10;
wp=1;
w1=1;
w2=1;
wz1=1;
wz2=1;
wd=1;
%init described in the paper as "Calibrated reconstruction"
%F_init=uu2F({u,v},{'None' '|F|=0'});
F_init=F;
F=correctF(F_init,f_prior,p_prior);
newF=reshape(F,1,9);
p=p_prior;
f=f_prior;
% options for lsqnonlin
o=optimoptions('lsqnonlin','Display','off');
%o=optimoptions('lsqnonlin','SpecifyObjectiveGradient',true,'Display','off');
o.Algorithm = 'levenberg-marquardt';
x=[newF p'];
[x,resnorm,residual,exitflag,output]=lsqnonlin(@(x) optimizex(x,u,v,wF^2,p_prior',wp,f_prior.^2,f_min.^2,f_ratio^2,w1,w2,wd,wz1,wz2), x.^2,[],[],o);
output.exitflag=exitflag;
F=reshape(x(1:9),3,3);
p=x(10:13);
p={p(1:2)' p(3:4)'};
f=F2f1f2(F,p);
end

function err=optimizex(x,u,v,wF,p_prior,wp,f_prior,f_min,f_ratio,w1,w2,wd,wz1,wz2)
F=reshape(x(1:9),3,3);
p=x(10:13);
f=F2f1f2(F,{p(1:2)' p(3:4)'}).^2;
errf=[w1*(f(1)-f_prior(1)) w2*(f(2)-f_prior(2)) wd*(f(1)*f_ratio-f(2))  wz1*(f(1)-f_min)*(f(1)<f_min)  wz2*(f(2)-f_min)*(f(2)<f_min)];
errp=wp*(p-p_prior);
errF=wF*sampson(F,u,v);
size(errf);
size(errp);
size(errF);
err=[errF errp errf];
% compute derivation of f?
end




