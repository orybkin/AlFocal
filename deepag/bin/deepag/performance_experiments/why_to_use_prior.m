% The script should not have existed. The theory behind it is plain wrong.
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

function [recomputed,original,f]=why_to_use_prior(file,corr,pop_size)
if nargin < 1
    global lots_of_inliers % determines whether we supposee we will use the methods in ransac or after we found a set of inliers
    lots_of_inliers=true;
    noise=1;
    file=['../../data/correspondences_F_synth_1K_' num2str(noise) 'noise.mat'];
    corr=40;
    pop_size=100;
    sizes=[0.001 0.01 0.1 0.2 0.5 1 2 5 10 100 1000]';
    priors=repmat([1500 2000],size(sizes,1),1).*sizes;
    colors={'r' 'g' 'k' 'c' '--r' '--g' '--k' '--c'  '-.r' '-.g' '-.k' '-.c'};
end

[recomputed,original,f]=calcFocals(file,corr,pop_size,priors);

img_idx=(abs(imag(f(:,1)))>eps) & (abs(imag(f(:,2)))>eps);

% histogram errors
figure();
hold on
cumhist(sort(original),20,1,'b');
for j=1:size(sizes,1)
    cumhist(sort(recomputed(:,j)),20,1,colors{j});
end
xlabel('error');
ylabel('frequency')
title('Given an estimate f, find best E. Numbers in legend are how off is the given f')
sizes_str=arrayfun(@(x) {['*' num2str(x)]},sizes);
legend([{'reference'}; sizes_str]) % char(176)]
hold off

end

function [recomputed,original,f]=calcFocals(file,corr,n, priors)
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
n_priors=size(priors,1);

recomputed=nan(n,n_priors);
original=nan(n,1);
support=zeros(n,1);
rng(867954152); 
global lots_of_inliers
for i=1:n
      repS = adprintf({}, [num2str(i), '/', num2str(n)]);
    %reshape
    uvector=u(:,i);
    points=size(uvector,1)/4;
    u1=reshape(uvector(1:end/2), 2, points);
    u2=reshape(uvector(end/2+1:end), 2, points);
    if ~lots_of_inliers
        %truncate
        sample=1:corr;
        sample=randperm(size(u1,2),corr);
        testsample=setdiff(1:size(u1,2),sample);
        testset={u1(:,testsample) u2(:,testsample)};
        u1=u1(:,sample);
        u2=u2(:,sample);
    else
        testset={nan};
    end
    % calculate estimation
    [F,A,support(i)]=F_features(u1,u2,'|F|=0',testset,3);
    F=reshape(F,3,3);
    f(i,:)=F2f1f2(F);
    f(i,:)=f(i,:)*diag([1/A{1}(1) 1/A{2}(1)]).*norm_(i,:);
    
    for j=1:n_priors
        % recomputed
        K1=diag([priors(j,1) priors(j,1) 1]);
        K2=diag([priors(j,2) priors(j,2) 1]);
        u1t=inv(K1)*a2h(u1);
        u2t=inv(K2)*a2h(u2);
        %testt={inv(K1)*a2h(testset{1}) inv(K2)*a2h(testset{2})};
        E_rec=peig5pt(u1t,u2t);
        E_rec=cell2mat(reshape(E_rec,1,1,[]));%{[testset{1} u1] [testset{2} u2]}
        E_rec=inliers2F(E_rec,nan,{u1t u2t});
        F_rec=inv(K2')*E_rec*inv(K1);

        % error measure
        % change pixels

        if lots_of_inliers
            test1=u1;
            test2=u2;
        else
            test1=testset{1};
            test2=testset{1};
        end
        recomputed(i,j)=sampson_error(F_rec,test1,test2);
    end
    
    original(i)=sampson_error(F,test1,test2);
  
    
    truth(i,:)=truth(i,:).*norm_(i,:);
    rmprintf(repS);
end

end