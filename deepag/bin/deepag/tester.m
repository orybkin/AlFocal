sceneType = {'randombox' 'random'};
pixel = 1/1000;
noise = 1*pixel; 
Npoints = 20;
Ncams = 2;
f1 = 3000*pixel;
f2 = 4000*pixel;
p1= [20*pixel 10*pixel];
p2= [20*pixel 10*pixel];
fmax=max(f1,f2);
Kgt1 = diag([f1 f1 1]);
Kgt2 = diag([f2 f2 1]);
Kgt1(1:2,3)=p1;
Kgt2(1:2,3)=p2;
gtk1 = 0;
gtk2 = 0;
[Pgt M m mgt] = GenerateScene(Npoints, [10*fmax 10*fmax], Ncams, 20*fmax, 30*fmax, 0, noise, [Kgt1;Kgt2], sceneType, [], [], [gtk1,gtk2], true);
Kgt{1} = Kgt1; Kgt{2}=Kgt2;
%ShowCameras(Pgt, Kgt, m, M, true, false, true, 1:7, mgt);
sample=randperm(size(m{1},2),7);
testsample=setdiff(1:Npoints,sample);
u1=m{1}(:,sample);
u2=m{2}(:,sample);
testset={m{1}(:,testsample) m{2}(:,testsample)};
[F,A]=F_features(u1,u2,'|F|=0',testset,0.001);
F=reshape(F,3,3);
f=F2f1f2(F)
F=F*(1+rand(3,3)/10);


[F,A]=uu2F_sampson(u1,u2);
f=F2f1f2(F)