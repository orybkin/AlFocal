function err=fcam5_experiment()
% Scene
X = [0 1 1 0 0 1 2 0 1 2
     0 0 1 1 0 0 1 1 1 2
     0 0 0 0 1 2 1 2 1 2];
X = X + 0.1*randn(size(X));
X = [X 2*rand(3,40)+[0;0;0]*ones(1,40)];
% Cameras

for i=1:5
    p{i}=[500; 500];
end
[K,R,C]=get_cameras(p);

err=fcam5_reconstruct(K,R,C,X);
err=norm(err);
end

function  [K,R,C]=get_cameras(p0)
% f
f0{1} = 1000;
f0{2} = 2000;
f0{3} = 3000;
f0{4} = 4000;
f0{5} = 2000;
% 1
R{1} = a2r([1;10;0],20/180*pi);
C{1} = 3*[2;1;-2];
% 2 
R{2} = a2r([7;10;2],-20/180*pi);
C{2} = 3*[-1;1;-2];
% 3
R{3} = a2r([10;5;1],30/180*pi);
C{3} = 3*[1;-1;-2];
if true   
    % 4
    R{4} = a2r([10;-10;1],-10/180*pi);
    C{4} = 5*[0.5;0.5;-2];
    % 5
    R{5} = a2r([-3;45;-3],-75/180*pi);
    C{5} = 3*[-2;0;0];
end
% K
for i=1:size(C,2)
    K{i} = [f0{i} 0 p0{i}(1);0 f0{i} p0{i}(2);0 0 1];
end

end
