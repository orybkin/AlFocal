-- Relative camera poses with two unknown focal length
-- x2^T K2^-T E K1^-1 x1 = 0, unknown f1, f2
--            F = K2^-T E K1^-1 x1
--            E = K2^T F K1 
-- T. Pajdla (pajdla@cvut.cz)

---- From the algebraic point of view, there seem not to be a difference 
---- between
---- K = [f 0 0; 0 f 0; 0 0 1] and K = [k11 k12 k13; 0 k22 k23; 0 0 1]
---- as the following code shows that there are no more algebraic 
---- constraints on F after eliminating f1, f2 than after eliminating 
---- k11, k12, k13, k22, k23,

---- K1, K2 diagonal
clearAll
R = QQ[f1,f2,f11,f12,f13,f21,f22,f23,f31,f32,f33,t1,t2,t3,t4,t5,t6,t7,t8,t9]
F = matrix{{f11,f12,f13},{f21,f22,f23},{f31,f32,f33}}
K1 = matrix{{f1, 0, 0}, {0, f1, 0}, {0, 0, 1}}
K2 = matrix{{f2, 0, 0}, {0, f2, 0}, {0, 0, 1}}
E = transpose(K2)*F*K1 -- Essential matrix
G = ideal(det(E)) + minors(1, 2*E*transpose(E)*E - trace(E*transpose(E))*E); -- Epipolar condition
dim G, codim G, degree G
H = flatten entries gens minors(2,F) -- rank of F must be 2
t = {t1,t2,t3,t4,t5,t6,t7,t8,t9};
H = apply({0,1,2,3,4,5,6,7,8},i->t#i*H#i)
H = ideal(fold(plus,0,H)+1)
dim H, codim H, degree H
Fundamental = eliminate(G,t)
Fundamental == G
Gs = saturate(G+H,ideal(f33*det(K1)*det(K2))); -- det(K1),det(K2) are non-zero
=ydim Gs, codim Gs, degree Gs
6M = eliminate(Gs,{f1,f2}); -- Constraints on F only
dim M, codim M, degree M
m = mingens M -- = -det(F)
g = mingens gb G;
s2 = mingens gb eliminate(Gs,f1)
s1 = mingens gb eliminate(Gs,f2)
-- Formula for f1
(m11,c11) = coefficients(s1_1_0,Variables=>{f1}) -- extract coefficients 
(m12,c12) = coefficients(s1_2_0,Variables=>{f1}) -- extract coefficients
(m13,c13) = coefficients(s1_3_0,Variables=>{f1}) -- extract coefficients
-- Formula for f2
(m21,c21) = coefficients(s2_1_0,Variables=>{f2}) -- extract coefficients 
(m22,c22) = coefficients(s2_2_0,Variables=>{f2}) -- extract coefficients
(m23,c23) = coefficients(s2_3_0,Variables=>{f2}) -- extract coefficients
-- Bougnoux formula for f1, f2
xx = x -> matrix{{0,- x_2,x_1},{x_2,0,- x_0},{- x_1,x_0,0}}
p1 = matrix{{0},{0},{1}}
p2 = matrix{{0},{0},{1}}
e1 = vector({f12*f23-f13*f22,-f11*f23+f13*f21,f11*f22-f12*f21})
e2 = vector({f21*f32-f31*f22,-f11*f32+f31*f12,f11*f22-f21*f12})
II = matrix({{1,0,0},{0,1,0},{0,0,0}})
nm1 = (-transpose(p1)*xx(e1)*II*transpose(F)*(p2*transpose(p2))*F*p1)_0_0
dn1 = (transpose(p1)*xx(e1)*II*transpose(F)*II*F*p1)_0_0
--- The 4th equation in Gs gives the Bougnoux Formula for first two rows of F independent
c23_0_1 == -nm1
c23_0_0 == dn1



---- K1 = I, K2 diagonal
clearAll
R = QQ[f2,f11,f12,f13,f21,f22,f23,f31,f32,f33, MonomialOrder=>Lex]
F = matrix{{f11,f12,f13},{f21,f22,f23},{f31,f32,f33}}
K2 = matrix{{f2, 0, 0}, {0, f2, 0}, {0, 0, 1}}
E = transpose(K2)*F -- Essential matrix
G = ideal(det(E)) + minors(1, 2*E*transpose(E)*E - trace(E*transpose(E))*E); -- Epipolar condition
dim G, codim G, degree G
Gs = saturate(G,ideal(det(K2))); -- det(K1) is non-zero
dim Gs, codim Gs, degree Gs
M = eliminate(Gs,{f2}); -- Constraints on F only
dim M, codim M, degree M
m = mingens M -- = -det(F)
gs = mingens gb Gs;
-- Formulae for f1
cofs = f->coefficients(f,Variables=>{f2})
gse = flatten entries gs;
cofsg = apply(gse,cofs);
cofsg_4
cofsg_5
cofsg_6
cofsg_7 -- this is formula (28)
cofsg_8
cofsg_9
cofsg_10
cofsg_11
cofsg_12

---- K1 = I, K2 diagonal, E = E^T = camera motion os a pure translation
clearAll
R = QQ[f2,f11,f12,f13,f21,f22,f23,f31,f32,f33, MonomialOrder=>Lex]
F = matrix{{f11,f12,f13},{f21,f22,f23},{f31,f32,f33}}
K2 = matrix{{f2, 0, 0}, {0, f2, 0}, {0, 0, 1}}
E = transpose(K2)*F -- Essential matrix
G = ideal(det(E)) + minors(1, 2*E*transpose(E)*E - trace(E*transpose(E))*E); 
G = ideal(minors(E+transpose(E))); -- Epipolar condition
dim G, codim G, degree G
Gs = saturate(G,ideal(det(K2))); -- det(K1) is non-zero
dim Gs, codim Gs, degree Gs
M = eliminate(Gs,{f2}); -- Constraints on F only
dim M, codim M, degree M
m = mingens M; -- = -det(F)
gs = mingens gb Gs;
-- Formulae for f1
cofs = f->coefficients(f,Variables=>{f2})
gse = flatten entries gs;
cofsg = apply(gse,cofs);
cofsg_1
cofsg_2
cofsg_16
cofsg_17
cofsg_19
cofsg_20
cofsg_21
cofsg_22
cofsg_23
cofsg_24

---- K1, K2 upper triangular
clearAll
R = QQ[a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,f11,f12,f13,f21,f22,f23,f31,f32,f33]
F = matrix{{f11,f12,f13},{f21,f22,f23},{f31,f32,f33}}
K1 = matrix{{a1, a2, a3}, {0, a4, a5}, {0, 0, 1}}
K2 = matrix{{b1, b2, b3}, {0, b4, b5}, {0, 0, 1}}
E = transpose(K2)*F*K1 -- Essential matrix
G = ideal(det(E)) + minors(1, 2*E*transpose(E)*E - trace(E*transpose(E))*E); -- Epipolar condition
Gs = saturate(G,ideal(det(K1)*det(K2))); -- det(K1),det(K2) are non-zero (takes apx 20 s)
M = eliminate(Gs,{a1,a2,a3,a4,a5,b1,b2,b3,b4,b5}); -- Constraints on F only
dim M, codim M, degree M
m = mingens M -- = -det(F)
---- We see that we have derived the same constraint on F, i.e. det(F)=0
---- However, these two cases are definitely different in the number of real 
---- solutions. Thus, we should see different discriminants of those two systems.

-- Parameterization with f->w is simpler
clearAll
R = QQ[f11,f12,f13,f21,f22,f23,f31,f32,f33,w1,w2, MonomialOrder=>Lex]
F = matrix{{f11,f12,f13},{f21,f22,f23},{f31,f32,f33}}
W1 = matrix{{1, 0, 0}, {0, 1, 0}, {0, 0, w1}}
W2 = matrix{{1, 0, 0}, {0, 1, 0}, {0, 0, w2}}
Ew = transpose(W2)*F*W1 -- Epipolar condition (K2/f2)^T*E*(K1/f1) and w1=1/f1, w2=1/f2 excludes f1=0 and f2=0
G =  ideal(det(Ew)) + minors(1, 2*Ew*transpose(Ew)*Ew - trace(Ew*transpose(Ew))*Ew);
Gs = saturate(G,ideal(w1*w2));
M = eliminate(Gs,{w1,w2});
dim M, codim M, degree M
m = mingens M

--- Two-sided Bougnoux formula for K1=K2 diagonal
clearAll;
R = QQ[f,f11,f12,f13,f21,f22,f23,f31,f32,f33] --- , MonomialOrder=>Lex]
F = matrix{{f11,f12,f13},{f21,f22,f23},{f31,f32,f33}};
K = matrix{{f, 0, 0}, {0, f, 0}, {0, 0, 1}};
E = K*F*K -- Essential matrix
G = ideal(det(E)) + minors(1, 2*E*transpose(E)*E - trace(E*transpose(E))*E); -- Epipolar condition
Gs = saturate(G,ideal(f)); -- det(K1)=det(K2) are non-zero
gs = mingens gb Gs;
-- Formulae for f
cofs = g->coefficients(g,Variables=>{f});
gse = flatten entries gs;
cofsg = apply(gse,cofs);
cofsg_2
--- a more conscise version
clearAll;
R = QQ[f,f11,f12,f13,f21,f22,f23,f31,f32,f33]
F = matrix{{f11,f12,f13},{f21,f22,f23},{f31,f32,f33}};
K = matrix{{f, 0, 0}, {0, f, 0}, {0, 0, 1}};
E = K*F*K;
G = ideal(det(E))+minors(1,2*E*transpose(E)*E-trace(E*transpose(E))*E);
Gs = saturate(G,ideal(f));
gse = flatten entries mingens gb Gs;
cofs = g->coefficients(g,Variables=>{f});
cofsg = apply(gse,cofs);
cofsg_2
--- the same by Joe
restart
R = QQ[e11,e12,e13,e21,e22,e23,e31,e32,e33,w,x11,x12,x13,x21,x22,x23,x31,x32,x33]
E = matrix{{e11,e12,e13},{e21,e22,e23},{e31,e32,e33}}
Kinv = matrix{{w,0,0},{0,w,0},{0,0,1}}
F = matrix{{x11,x12,x13},{x21,x22,x23},{x31,x32,x33}}
I = minors(1,F - Kinv*E*Kinv) + minors(1,2*E*transpose(E)*E - trace(E*transpose(E))*E) + ideal(det(E))
S = flatten entries gens eliminate({e11,e12,e13,e21,e22,e23,e31,e32,e33},I)
toString S#2
factor S#2


