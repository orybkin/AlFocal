% The script shows how to use absolute value of the focal length
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

function [best,estion,support, baseline,truth]=how_to_use_abs(file,corr,pop_size, method, noise_out)
if nargin < 1
    file='../../data/paris/correspondences_F_synth_1K_1noise.mat';
    corr=7;
    pop_size=100;
    method='Prop6';
    noise_out=10;
end

[corrected,recomputed,reransaced]=calcFocals(file,corr,pop_size);

basedata=get_foc_error(baseline,truth);
estdata=get_foc_error(estion,truth);
bestdata=get_foc_error(best,truth);

% histogram errors
figure();
cumhist(sort(estdata),20,1,'-g');
hold on
cumhist(sort(basedata),20,1,'-r');
cumhist(sort(bestdata),20,1,'-b');
hold off
legend('6-Ratio','7pt','best among 6-Ratio')
end

function [corrected,recomputed,reransaced]=calcFocals(file,corr,n)
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
f=zeros(n,2);
corrected=nan(n);
recomputed=nan(n);
reransaced=nan(n);
support=zeros(n,1);
rng(867954152); 


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
    
    % calculate estimation
    [F,A,support(i)]=F_features(u1,u2,'|F|=0',testset,3,false);
    f(i,:)=F2f1f2(reshape(F,3,3));
    f(i,:)=f(i,:)*diag([1/A{1}(1) 1/A{2}(1)]).*norm_(i,:);
    
    % corrected
    corrF=correctF(F,f(i,:));
    corrected(i)=norm(sampson(coorF,u1,u2));
    
    % recomputed
    TODO
    
    % reransaced (optional)    
    
    truth(i,:)=truth(i,:).*norm_(i,:);
    rmprintf(repS);
end
end