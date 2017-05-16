% The script histograms Prior solver results versus baseline
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

function [hartley_ratio,hartley,prior,baseline]=hartley_ratio_perf(corr,pop_size, noise)
if nargin < 1
    noise=1;
    corr=7;
    pop_size=100;
end

[hartley_ratio,hartley,prior,baseline]=calcFocals(noise,corr,pop_size);

% histogram iterations
figure();
hold on
histhist(sort(hartley_ratio.iter),15,1,'-g');
histhist(sort(hartley.iter),15,1,'-r');
hold off
xlabel('iterations')
ylabel('frequency')
legend('hartley ratio iterations', 'hartley iterations')
saveas(gcf,'hartley_ratio_iterations.fig')

% histogram f1
figure();
hold on
histhist(sort(hartley_ratio.f(:,1)),15,1,'-g');
histhist(sort(hartley.f(:,1)),15,1,'-r');
hold off
xlabel('f1')
ylabel('frequency')
legend('hartley ratio f', 'hartley f')
saveas(gcf,'hartley_ratio_f1.fig')


% histogram f2
figure();
hold on
histhist(sort(hartley_ratio.f(:,2)),15,1,'-g');
histhist(sort(hartley.f(:,2)),15,1,'-r');
hold off
xlabel('f2')
ylabel('frequency')
legend('hartley ratio f', 'hartley f')
saveas(gcf,'hartley_ratio_f2.fig')

truth=repmat([3 4],size(hartley_ratio.f,1),1);
ratio_error=get_foc_error(hartley_ratio.f,truth);
hartley_error=get_foc_error(hartley.f,truth);
% histogram errors
figure();
hold on
cumhist(sort(ratio_error),15,'-g');
cumhist(sort(hartley_error),15,'-r');
hold off
legend('hartley ratio error', 'hartley error')

%toc();
end

function [hartley_ratio,hartley,prior,baseline]=calcFocals(noise,corr,n,prior )
global debugg;
%get n focal length estimations of the bougnoux formula on corr coordinates from the data in
%specified file
prior.f=zeros(n,2);
prior.p1=zeros(n,2);
prior.p1=zeros(n,2);

hartley_ratio.f=zeros(n,2);
hartley_ratio.p1=zeros(n,2);
hartley_ratio.p2=zeros(n,2);
hartley_ratio.iter=zeros(n,1);
hartley_ratio.converge=zeros(n,1);

hartley.f=zeros(n,2);
hartley.p1=zeros(n,2);
hartley.p2=zeros(n,2);
hartley.iter=zeros(n,1);
hartley.converge=zeros(n,1);

baseline=zeros(n,2);
rng(867954152); 


for i=1:n
    repS = adprintf({}, [num2str(i), '/', num2str(n)]);
    %generate
    sceneType = {'randombox' 'random'};
    pixel = 1/1000;
    noise = noise*pixel; 
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
    sample=randperm(size(m{1},2),corr);
    testsample=setdiff(1:Npoints,sample);
    u1=m{1}(:,sample);
    u2=m{2}(:,sample);
    testset={m{1}(:,testsample) m{2}(:,testsample)};
    
    
    %calculate baseline estimation
    [F,A]=F_features(u1,u2,'|F|=0',testset,0.001);
    F=reshape(F,3,3);
    baseline(i,:)=F2f1f2(F);
    ratio=abs(F2ratio(F));
    real_idx=@(x) ((imag(x(:,1))<eps) & (imag(x(:,2)))<eps);
    img_idx=@(x) ((imag(x(:,1))>eps) & (imag(x(:,2)))>eps);
    if ~real_idx(baseline(i,:)) && ~img_idx(baseline(i,:))
        ratio=1;
    end
    
    baseline(i,:)=baseline(i,:)*diag([1/A{1}(1) 1/A{2}(1)]);
    p1=[0;0];
    p2=[0;0];
    
    % Priors
    f1prior =f1+ 1.5*rand(1);
    u1prior = p1(1)+0.1*rand(1);
    v1prior = p1(2)+0.1*rand(1);
    
    f2prior = f2+ 1.5*rand(1);
    u2prior = p2(1)+0.05*rand(1);
    v2prior = p2(2)+0.05*rand(1);
    prior.f(i,:)=[f1prior f2prior];
    prior.p1(i,:)=[u1prior v1prior];
    prior.p2(i,:)=[u2prior v2prior];
    % Hartley with ratio
    tic()
    [~, p_hartley, f_hartley, output] = uu2F_hartley(F,u1,u2,[f1prior; f2prior],{[u1prior; v1prior] [u1prior; v1prior]}, 0, ratio);
    toc()
    hartley_ratio.converge(i)=output.exitflag;
    hartley_ratio.iter(i)=output.iterations;
    hartley_ratio.f(i,:)=f_hartley;
    hartley_ratio.p1(i,:)=p_hartley{1}';
    hartley_ratio.p2(i,:)=p_hartley{2}';
    norm(f_hartley-[f1 f2])
    
    % Hartley
    tic()
    [~, p_hartley, f_hartley, output] = uu2F_hartley(F,u1,u2,[f1prior; f2prior],{[u1prior; v1prior] [u1prior; v1prior]});
    toc()
    hartley.converge(i)=output.exitflag;
    hartley.iter(i)=output.iterations;
    hartley.f(i,:)=f_hartley;
    hartley.p1(i,:)=p_hartley{1}';
    hartley.p2(i,:)=p_hartley{2}';
    norm(f_hartley-[f1 f2])
    
    rmprintf(repS);
    end

end