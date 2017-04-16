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
R = QQ[f1,f2,f11,f12,f13,f21,f22,f23,f31,f32,f33, MonomialOrder=>Lex]
F = matrix{{f11,f12,f13},{f21,f22,f23},{f31,f32,f33}}
K1 = matrix{{f1, 0, 0}, {0, f1, 0}, {0, 0, 1}}
K2 = matrix{{f2, 0, 0}, {0, f2, 0}, {0, 0, 1}}
E = transpose(K2)*F*K1 -- Essential matrix
IF = ideal det F
dim IF, codim IF, degree IF
G = ideal(det(E)) + minors(1, 2*E*transpose(E)*E - trace(E*transpose(E))*E); -- Epipolar condition
dim G, codim G, degree G
Gs = saturate(G,ideal(f1*f2)); -- det(K1),det(K2) are non-zero
-- Gs = Gs + ideal(f33-1)
dim Gs, codim Gs, degree Gs
M = eliminate(Gs,{f1,f2}); -- Constraints on F only
dim M, codim M, degree M
m = mingens M -- = -det(F)
g = mingens gb Gs
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

nm2 = (-transpose(p2)*xx(e2)*II*F*(p1*transpose(p1))*transpose(F)*p2)_0_0
dn2 = (transpose(p2)*xx(e2)*II*F*II*transpose(F)*p2)_0_0

factor c11_0_1*c12_0_0-c11_0_0*c12_0_1 -- not equal!

fl1=1500
fl2=2000
K1 = matrix{{fl1,0,0},{0,fl1,0},{0,0,1}}
K2 = matrix{{fl2,0,0},{0,fl2,0},{0,0,1}}
R2 = matrix{{0,1,0},{-3/5,0,4/5},{4/5,0,3/5}} * matrix{{1,0,0},{0,3/5,-4/5},{0,4/5,3/5}} ** R
R1 = id_(R^3)
points=7
X  = matrix(fillMatrix(mutableMatrix(QQ,4,points)))
X =  mutableMatrix for i from 0 to points-1 list (1/X_i_3)*(X_i)
for i from 0 to points-1 do X_(2,i)=X_(2,i)+3
X = matrix X
t1 = transpose matrix{{0,0,0}}
t2 = transpose matrix{{-1,-1,-1}}
P1=K1*(R1| -R1*t1)
P2=K2*(R2| -R2*t2)

u1=P1*X
u2=P2*X

A = transpose matrix for i from 0 to points-1 list u1_i ** u2_i 

--  A = matrix(fillMatrix(mutableMatrix(ZZ,8,9)))
m = transpose matrix {{f11,f12,f13,f21,f22,f23,f31,f32,f33}}
rank A
eA = A*m
k = kernel A
Fm=matrix{{k_0_0,k_0_1,k_0_2},{k_0_3,k_0_4,k_0_5},{k_0_6,k_0_7,k_0_8}}
IA = minors(1,eA)
GAs = Gs + IA
gGAS=mingens gb GAs
dim GAs, codim GAs, degree GAs
GAs = GAs+ideal(f33-1)
pGAS = minimalPrimes GAs

-- degenerace? pohyb nekam
-- Deg= Gs + ideal(f11,f12,f13,f21,f22,f23,f31,f32,f33)
Deg3=ideal(f11,f12+1,f13,f21-1,f22,f23,f31,f32,f33)
mingens gb De



GA = G + IA
mingens gb GA
dim GA, codim GA, degree GA
