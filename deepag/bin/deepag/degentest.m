% The script gives an example of degeneration in the sense of 
% [On the existence of Epipolar Matrices], theorem 2.2
% turns out the the situation described there for 7 points does not define
% 7 linearly independent constraints - only 6
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

points=7
% focal lengths
f1 = 900;
f2 = 1100;
sc = max(f1,f2);
% 3D points
points=7;
X=randn(3,points);
dir=randn(3,1); % projects 6 points to plane
dir=dir/norm(dir);
X(:,1:end-1)=(eye(3)-dir*dir')*X(:,1:end-1);
% Internal camera calibration
K1 = [f1  0    0
    0   f1   0
    0   0    1];
K2 = [f2 0  0
    0  f2 0
    0  0  1];
% Rotations
R1 = a2r([1 0 0]',-pi/10);
R2 = a2r([1 0 0]',pi/11);
% Projection matrices
P1 = K1*R1*[eye(3)     [0;0;0]];
P2 = K2*R2*[eye(3)  [0;sc;0]];
% Image points
x{1} = X2u(X,P1);
x{2} = X2u(X,P2);


% x{1}=randn(2,7);
% x{2}=[[linspace(randn(1,1), randn(1,1),6);
% linspace(randn(1,1), randn(1,1),6)] rand(2,1)];
% x{1}=a2h(x{1});
% x{2}=a2h(x{2});
[F,A]=F_features(x{1},x{2},'|F|=0')
Ft=reshape(F,3,3);
if true
figure()
scatter(x{2}(1,:),x{2}(2,:),1)
figure()
scatter(x{1}(1,:),x{1}(2,:),1)
end
B = zeros(size(x{1},2),9);
% B * f = 0
for i=1:size(x{1},2)
    B(i,:) = m2v(x{1}(:,i)*x{2}(:,i)');
end