function F=correctF(F,f,p)
% Corrects F to be consistent with f and p, and still be a fundamental
% matrix. Described in the paper
if nargin<3
    p = {[0;0] [0;0]};
end
K1=diag([f(1) f(1) 1]);
K2=diag([f(2) f(2) 1]);
K1(1:2,3)=p(1:2);
K2(1:2,3)=p(3:4);

E=K2'*F*K1;
[U,D,V]=svd(E);
E=U*diag([1 1 0])*V'; % an interesting point is that this step is not required. maybe try to leave it out and test converge
F=inv(K2')*E*inv(K1);
end
