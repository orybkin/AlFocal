A = rand(3,3);

[U, D, V] = svd(A);
R = U';
T = [0 3 -2; -3 0 4; 2 -4 0];

f1gt = 3;
u1gt = 0.01;
v1gt = 0.02;

f2gt = 4;
u2gt = 0.02;
v2gt = 0.001;

K1 = [f1gt, 0 u1gt; 0 f1gt v1gt; 0 0 1];
K2 = [f2gt, 0 u2gt; 0 f2gt v2gt; 0 0 1];
F1 = transpose(K1^(-1))*T*R*(K1^(-1));
F2 = transpose(K2^(-1))*T*R*(K1^(-1));

f1prior =f1gt+ 1.5*rand(1);
u1prior = u1gt+0.01*rand(1);
v1prior = v1gt+0.01*rand(1);
w1 =0.1;
w2 =1;

f2prior = f2gt+ 1.5*rand(1);
u2prior = u2gt+0.01*rand(1);
v2prior = v2gt+0.01*rand(1);
w3 =0.1;
w4 =1;

if 0
    % one focal
    [fo, uo, vo, l1, l2, err, iter] = focal_from_F(F1, f1prior, u1prior, v1prior, w1, w2)

    Ko = [fo 0 uo; 0 fo vo; 0 0 1]
    Eo = transpose(Ko)*F1*Ko;
    [Uo Do Vo] = svd(Eo);
    Do
end

[f1o, u1o, v1o, f2o, u2o, v2o, l1, l2, err, iter] = f1_f2_from_F(F2, f1prior, u1prior, v1prior,f2prior, u2prior, v2prior, w1, w2, w3, w4);

K1o = [f1o 0 u1o; 0 f1o v1o; 0 0 1];
K2o = [f2o 0 u2o; 0 f2o v2o; 0 0 1];

E2o = transpose(K2o)*F2*K1o;
[U2o D2o V2o] = svd(E2o);
D2o;


