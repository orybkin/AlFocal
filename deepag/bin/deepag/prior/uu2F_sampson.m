function [ F, output ] = uu2F_sampson(u, v)
% UU2F_HARTLEY solves 7pt problem as described in paper [Reconstruction
% from two views using approximate calibration]
% transformation
if size(u,1)<3
    u=[u; ones(1,size(u,2))];
    v=[v; ones(1,size(u,2))];
end
%init described in the paper as "Calibrated reconstruction"
[F_init,~]=F_features(u,v,'|F|=0');
F=reshape(F_init,1,9);
% options for lsqnonlin
o=optimoptions('lsqnonlin','Display','off');
o.Algorithm = 'levenberg-marquardt';
[F,resnorm,residual,exitflag,output]=lsqnonlin(@(F) optimizeF(F,u,v), F,[],[],o);
output.exitflag=exitflag;
F=reshape(F(1:9),3,3);
end

function err=optimizeF(F,u,v)
F=reshape(F,3,3);
err=sampson(F,u,v);
% compute derivation of f?
end




