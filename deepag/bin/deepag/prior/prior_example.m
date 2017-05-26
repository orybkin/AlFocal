% comparison of Hartley v. Zuzana prior methods - an example 
%


%% generating
sceneType = {'randombox' 'random'};
pixel = 1/1000;
noise = 1*pixel;  %noise=0;
Npoints = 20;
Ncams = 2;
samplesize = 7;
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
sample=randperm(size(m{1},2),samplesize);
testsample=setdiff(1:Npoints,sample);
u={m{1}(:,sample) m{2}(:,sample)};
testset={m{1}(:,testsample) m{2}(:,testsample)};
%% calculating
% 7pt
clear estion3 estion estion6;
[F,A]=F_features(u{1},u{2},'|F|=0',testset,0.001);
F_t = inv(A{2})'*reshape(F,3,3)*inv(A{1});
f_7pt=F2f1f2(F_t)*diag([1/A{1}(1) 1/A{2}(1)]);
f_7pt
norm(abs(f_7pt)-[f1 f2])
%% calculating
% Zuzana
f1prior =f1+ 1.5*rand(1);
u1prior = p1(1)+0.1*rand(1);
v1prior = p1(2)+0.1*rand(1);
w1 =0.1;
w2 =1;

f2prior = f2+ 1.5*rand(1);
u2prior = p2(1)+0.1*rand(1);
v2prior = p2(2)+0.1*rand(1);
w3 =0.1;
w4 =1;
tic()
[f1_zuz, u1o, v1o, f2_zuz, u2o, v2o, l1, l2, err, iter] = f1_f2_from_F(F_t, f1prior, u1prior, v1prior,f2prior, u2prior, v2prior, w1, w2, w3, w4);
toc()
norm([f1_zuz f2_zuz]-[f1 f2])

%% Hartley
tic()
[F, p_hartley, f_hartley, output] = uu2F_hartley(F_t,u{1},u{2},[f1prior; f2prior],{[u1prior; v1prior] [u1prior; v1prior]});
toc()
output.iterations
norm(f_hartley-[f1 f2])
f_hartley
%%
norm([f1 f2]-[f1prior f2prior])


