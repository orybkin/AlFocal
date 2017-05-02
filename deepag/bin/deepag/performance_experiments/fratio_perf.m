% The script histograms 

% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

function [pol1,pol,boug,truth]=fratio_perf(file,corr,pop_size, noise_out)
if nargin < 1
    noise=1;
    file=['../../data/paris/correspondences_F_synth_1K_' num2str(noise) 'noise.mat'];
    corr=7;
    pop_size=100;
    noise_out=0;
end

[pol1,pol,boug,truth]=calcFocals(file,corr,pop_size, noise_out);

pol1_e=get_rat_error_add(pol1,truth);
pol_e=get_rat_error_add(pol,truth);
boug_e=get_rat_error_add(boug,truth);

% histogram errors
figure();
cumhist(sort(pol1_e),20,1);
hold on
cumhist(sort(pol_e),20,1);
cumhist(sort(boug_e),20,1);
hold off
legend('f1/f2','direct','Bougnoux')

%toc();
end

function [pol1,pol,boug,truth]=calcFocals(file,corr,n,noise)
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
pol1=ones(n,2);
pol=ones(n,2);
boug=ones(n,2);
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
    
    pol1(i,1)=F2ratio(reshape(Fund,3,3),{[0;0] [0;0]},'Polynomial1');
    pol1(i,1)=pol1(i,1)*A{1}(1)/A{2}(1)*norm_(i,2)/norm_(i,1);
    
    pol(i,1)=F2ratio(reshape(Fund,3,3),{[0;0] [0;0]},'Polynomial');
    pol(i,1)=pol(i,1)*A{1}(1)/A{2}(1)*norm_(i,2)/norm_(i,1);
    
    boug(i,1)=F2ratio(reshape(Fund,3,3),{[0;0] [0;0]},'Bougnoux');
    boug(i,1)=boug(i,1)*A{1}(1)/A{2}(1)*norm_(i,2)/norm_(i,1);
    
    truth(i,:)=truth(i,:).*norm_(i,:);
    rmprintf(repS);
end
truth(:,1)=truth(:,2)./truth(:,1);
truth(:,2)=1;
end