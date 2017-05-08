% The script histograms Prior solver results versus baseline
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

function [zuzana,hartley,prior,baseline,truth]=hartley_perf(file,corr,pop_size, method, noise)
if nargin < 1
    noise=1;
    file=['../../data/correspondences_F_synth_1K_' num2str(noise) 'noise.mat'];
    corr=7;
    pop_size=100;
    method='Prop6';
    noise=1;
end

[zuzana,hartley,prior,baseline,truth]=calcFocals(file,corr,pop_size);

% histogram errors
figure();
hold on
cumhist(sort(zuzana.iter),20,1,'-g');
cumhist(sort(hartley.iter),20,1,'-r');
hold off
legend('zuzana iterations', 'hartley iterations')

%toc();
end

function [zuzana,hartley,prior,baseline,truth]=calcFocals(file,corr,n,prior )
global debugg;
%get n focal length estimations of the bougnoux formula on corr coordinates from the data in
%specified file

%load data
load(file);
u=[corr_tr.u corr_val.u corr_tst.u];
truth=[corr_tr.f corr_val.f corr_tst.f];
norm_=[corr_tr.norm corr_val.norm corr_tst.norm];
%truncate
u=u(:,1:n); truth=truth(:,1:n)'; norm_=norm_(:,1:n)';

prior.f=zeros(n,2);
prior.p1=zeros(n,2);
prior.p1=zeros(n,2);

zuzana.f=zeros(n,2);
zuzana.p1=zeros(n,2);
zuzana.p2=zeros(n,2);
zuzana.iter=zeros(n,1);
zuzana.converge=zeros(n,1);

hartley.f=zeros(n,2);
hartley.p1=zeros(n,2);
hartley.p2=zeros(n,2);
hartley.iter=zeros(n,1);
hartley.converge=zeros(n,1);

baseline=zeros(n,2);
rng(867954152); %the reason I create noise beforehand is that now I can be sure
%for every call of this method the noise will be the same, which allows for
%comparison


for i=1:n
    repS = adprintf({}, [num2str(i), '/', num2str(n)]);
    %reshape
    sceneType = {'randombox' 'random'};
    pixel = 1/1000;
    noise = 1*pixel;  noise=1;
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
    u1=m{1}(:,sample);
    u2=m{2}(:,sample);
    testset={m{1}(:,testsample) m{2}(:,testsample)};
    
    
    %calculate baseline estimation
    [F,A]=F_features(u1,u2,'|F|=0',testset,0.001);
    F=reshape(F,3,3);
    baseline(i,:)=F2f1f2(F);
    baseline(i,:)=baseline(i,:)*diag([1/A{1}(1) 1/A{2}(1)]).*norm_(i,:);
    f1=baseline(i,1);
    f2=baseline(i,2);
    p1=[0;0];
    p2=[0;0];
    
    % Zuzana
    f1prior =f1+ 1.5*rand(1);
    u1prior = p1(1)+0.1*rand(1);
    v1prior = p1(2)+0.1*rand(1);
    w1 = 0.1;
    w2 = 1;
    
    f2prior = f2+ 1.5*rand(1);
    u2prior = p2(1)+0.05*rand(1);
    v2prior = p2(2)+0.05*rand(1);
    prior.f(i,:)=[f1prior f2prior];
    prior.p1(i,:)=[u1prior v1prior];
    prior.p2(i,:)=[u2prior v2prior];
    w3 =0.1;
    w4 =1;
    tic()
    [f1_zuz, u1_zuz, v1_zuz, f2_zuz, u2_zuz, v2_zuz, l1, l2, err, iter] = f1_f2_from_F(F, f1prior, u1prior, v1prior,f2prior, u2prior, v2prior, w1, w2, w3, w4);
    toc()
    zuzana.converge(i)=iter>50;
    zuzana.iter(i)=iter;
    zuzana.f(i,:)=[f1_zuz f2_zuz];
    zuzana.p1(i,:)=[u1_zuz v1_zuz];
    zuzana.p2(i,:)=[u2_zuz v2_zuz];
    
    % Hartley
    tic()
    [F, p_hartley, f_hartley, output] = uu2F_hartley(F,u1,u2,[f1prior; f2prior],{[u1prior; v1prior] [u1prior; v1prior]});
    toc()
    hartley.converge(i)=output.exitflag;
    hartley.iter(i)=output.iterations;
    hartley.f(i,:)=f_hartley;
    hartley.p1(i,:)=p_hartley{1}';
    hartley.p2(i,:)=p_hartley{2}';
    norm(f_hartley-[f1 f2])
    truth(i,:)=truth(i,:).*norm_(i,:);
    
    rmprintf(repS);
end
end