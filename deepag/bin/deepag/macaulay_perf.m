% The script histograms 6-Ratio results versus baseline
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

function [pol1,pol2,pol3,boug,truth]=macaulay_perf(file,corr,pop_size, method, noise)
if nargin < 1
    file='../../data/paris/correspondences_F_synth_10K_1noise.mat';
    corr=7;
    pop_size=10000;
    method='Prop6';
    noise=1;
end

[pol1,pol2,pol3,boug,truth]=calcFocals(file,corr,pop_size, method, noise);

pol1=diffify(pol1,truth);
pol2=diffify(pol2,truth);
pol3=diffify(pol3,truth);
boug=diffify(boug,truth);

% histogram errors
figure();
cumhist(sort(pol1),20,1,'-g');
hold on
cumhist(sort(pol2),20,1,'-r');
cumhist(sort(pol3),20,1,'-b');
cumhist(sort(boug),20,1,'-k');
hold off
legend('1','2','3','Bougnoux')

%toc();
end

function estdata=diffify(estion,truth)
%transform data to difference of f-Ratio and ground truth f-Ratio (errors)

bigger=abs(estion(:,1))>abs(truth(:,1));
estdata(bigger)=abs(truth(bigger,1))./abs(estion(bigger,1));
estdata(not(bigger))=abs(estion(not(bigger),1))./abs(truth(not(bigger),1));
estdata=1-estdata;
end

function [pol1,pol2,pol3,baseline,truth]=calcFocals(file,corr,n,method,noise)
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
pol1=zeros(n,2);
pol2=zeros(n,2);
pol3=zeros(n,2);
baseline=zeros(n,2);
rng(867954152); %the reason I create noise beforehand is that now I can be sure
%for every call of this method the noise will be the same, which allows for
%comparison
noisy=noise*randn(2,n);
noisy(:,1);


for i=1:n
    repS = adprintf({}, [num2str(i), '/', num2str(n)]);
    %reshape
    uvector=u(:,i);
    points=size(uvector,1)/4;
    u1=reshape(uvector(1:end/2), 2, points);
    u2=reshape(uvector(end/2+1:end), 2, points);
    %truncate
    sample=1:corr;
    %sample=randperm(size(u1,2),corr);
    testsample=setdiff(1:size(u1,2),sample);
    testset={u1(:,testsample) u2(:,testsample)};
    u1=u1(:,sample);
    u2=u2(:,sample);
    %noise
    u1(:,1)=u1(:,1)+noisy(:,i);
    
    
    %calculate baseline estimation
    [Fund,A]=F_features(u1,u2,'|F|=0',testset,3,false);
    baseline(i,:)=F2f1f2(reshape(Fund,3,3));
    baseline(i,:)=baseline(i,:)*diag([1/A{1}(1) 1/A{2}(1)]).*norm_(i,:);
    
    pol1(i,:)=F2f1f2(reshape(Fund,3,3),{[0;0] [0;0]},'Polynomial1');
    pol1(i,:)=pol1(i,:)*diag([1/A{1}(1) 1/A{2}(1)]).*norm_(i,:);
    
    pol2(i,:)=F2f1f2(reshape(Fund,3,3),{[0;0] [0;0]},'Polynomial1');
    pol2(i,:)=pol2(i,:)*diag([1/A{1}(1) 1/A{2}(1)]).*norm_(i,:);
    
    pol3(i,:)=F2f1f2(reshape(Fund,3,3),{[0;0] [0;0]},'Polynomial3');
    pol3(i,:)=pol3(i,:)*diag([1/A{1}(1) 1/A{2}(1)]).*norm_(i,:);
    
    truth(i,:)=truth(i,:).*norm_(i,:);
    rmprintf(repS);
end
end