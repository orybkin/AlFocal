% The script histograms 6-Ratio results versus baseline
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

function [best,estion,support, baseline,truth]=ratio6_perf(file,corr,pop_size, method, noise_out)
if nargin < 1
    noise=1;
    file=['../../data/correspondences_F_synth_1K_' num2str(noise) 'noise.mat'];
    corr=7;
    pop_size=100;
    method='Prop6';
    noise_out=10;
end

[best, estion,support,baseline,truth]=calcFocals(file,corr,pop_size, method, noise_out);

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
    
if false
    figure();
%     histogram(abs(basedata),[1:30]/20);
%     hold on
%     plot(basedata,cumsum(sort(basedata)));
    %title('(quadratic) errors of 7pt');
    %xlabel('|f_{6-R}-f_{ground}');
    title('baseline f_1 histogram');
    xlabel('f_1');
    ylabel('frequency');
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperPosition', [0 0 1024 800]);
    
    figure();
%     histogram(abs(estdata),[1:30]/20);
%     hold on
%     plot(sort(estdata),cumsum(sort(estdata)));
    %title('(quadratic) errors of 6-Ratio');
    %xlabel('|f_{6-R}-f_{ground}');
    title('6-Ratio f_1 histogram');
    xlabel('f_1');
    ylabel('frequency');
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperPosition', [0 0 1024 800]);
end
%toc();
end

function [best,estion,support,baseline,truth]=calcFocals(file,corr,n,method,noise)
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
estion=zeros(n,2);
baseline=zeros(n,2);
best=zeros(n,2);
support=zeros(n,1);
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
    
    
    %calculate 6-Ratio estimation
    [Fund,A,support(i)]=F_features(u1,u2,method,testset,3,false);
    estion(i,:)=F2f1f2(reshape(Fund,3,3));
    estion(i,:)=estion(i,:)*diag([1/A{1}(1) 1/A{2}(1)]).*norm_(i,:);
    truth(i,:)=truth(i,:).*norm_(i,:);
    %calculate baseline estimation
    [Fund,A]=F_features(u1,u2,'|F|=0',testset,3,false);
    baseline(i,:)=F2f1f2(reshape(Fund,3,3));
    baseline(i,:)=baseline(i,:)*diag([1/A{1}(1) 1/A{2}(1)]).*norm_(i,:);
    truth(i,:)=truth(i,:).*norm_(i,:);
    %calculate best one
    [F6,A6,f6]  = uu2F({u1 u2},{'None','Prop6'},testset);
    
    choices=size(F6,3);
    estion6=zeros(choices,2);
    for j=1:choices
        F_t = inv(A{2})'*F6(:,:,j)*inv(A{1});
        estion6(j,:)=F2f1f2(F_t)*diag([1/A6{1}(1) 1/A6{2}(1)]);
        % differ=differ+norm(estion-estion6{j},'fro');
    end
    % what to do with imaginary things
    besti=abs(estion6)-repmat(truth(i,:),choices,1);
    besti=besti(:,1).^2+besti(:,2).^2;
    [~,besti]=min((besti));
    best(i,:)=estion6(besti,:);
    rmprintf(repS);
end
end