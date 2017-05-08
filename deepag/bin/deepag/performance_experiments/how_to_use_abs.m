% The script shows how to use absolute value of the focal length
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

function [corrected,recomputed,original,f]=how_to_use_abs(file,corr,pop_size, method, noise_out)
if nargin < 1
    global lots_of_inliers % determines whether we supposee we will use the methods in ransac or after we found a set of inliers
    lots_of_inliers=false;
    noise=1;
    file=['../../data/correspondences_F_synth_1K_' num2str(noise) 'noise.mat'];
    corr=7;
    pop_size=100;
    method='Prop6';
    noise_out=10;
end

[corrected,recomputed,original,weird,f]=calcFocals(file,corr,pop_size);

img_idx=(abs(imag(f(:,1)))>eps) | (abs(imag(f(:,2)))>eps);

% histogram errors
if false
    figure();
    hold on
    cumhist(sort(corrected(~img_idx)),20,1,'-g');
    cumhist(sort(recomputed(~img_idx)),20,1,'-r');
    cumhist(sort(original(~img_idx)),20,1,'-r');
    cumhist(sort(weird(~img_idx)),20,1,'-r');
    title('real')
    hold off
    legend('corrected','recomputed','original','E2PP')
end

figure();
hold on
cumhist(sort(original(~img_idx)),20,1,'-g');
cumhist(sort(corrected(img_idx)),20,1,'-g');
cumhist(sort(recomputed(img_idx)),20,1,'-r');
cumhist(sort(original(img_idx)),20,1,'-r');
cumhist(sort(weird(img_idx)),20,1,'-r')
xlabel('error');
ylabel('frequency')
title('How to use imaginary focal length')
hold off

legend('reference','corrected','recomputed','original','E2PP')
end

function [corrected,recomputed,original,weird,f]=calcFocals(file,corr,n)
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
corrected=nan(n,1);
recomputed=nan(n,1);
weird=nan(n,1);
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
        
    % corrected
    [F_corr,K1,K2]=correctF(F,abs(f(i,:)));
    if norm(F/F(1)-F_corr/F_corr(1),'fro')>1e-6 & (abs(imag(f(i,1)))<eps) & (abs(imag(f(i,2)))<eps)
        norm(F/F(1)-F_corr/F_corr(1),'fro')
        f(i,:)
        norm_(i,:)
        A{1}
        A{2}
        warning('calcFocals');
    end
    % recomputed
    u1t=inv(K1)*a2h(u1);
    u2t=inv(K2)*a2h(u2);
    %testt={inv(K1)*a2h(testset{1}) inv(K2)*a2h(testset{2})};
    E_rec=peig5pt(u1t,u2t);
    E_rec=cell2mat(reshape(E_rec,1,1,[]));%{[testset{1} u1] [testset{2} u2]}
    E_rec=inliers2F(E_rec,nan,{u1t u2t});
    F_rec=inv(K2')*E_rec*inv(K1);
    % method used by Tomas in frecexample
    PP = E2PP(K2'*F*K1);
    E_weird = PP2F(PP{1}{1},PP{2}{1});
    F_weird = inv(K2')*E_weird*inv(K1);
    
    % error measure
    % change pixels
    
    if lots_of_inliers
        test1=u1;
        test2=u2;
    else
        test1=testset{1};
        test2=testset{1};
    end
    original(i)=sampson_error(F,test1,test2);
    corrected(i)=sampson_error(F_corr,test1,test2);
    recomputed(i)=sampson_error(F_rec,test1,test2);    
    weird(i)=sampson_error(F_weird,test1,test2);
  
    % reransaced (optional)    
    
    truth(i,:)=truth(i,:).*norm_(i,:);
    rmprintf(repS);
end
end