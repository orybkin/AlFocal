% The script histograms 6-Ratio results versus baseline
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

function [percent7]=imaginary_number(file,corr,pop_size, method, noise)
if nargin < 1
    noise=1;
    file=['../../data/correspondences_F_synth_1K_' num2str(noise) 'noise.mat'];
    corrs=8:40;
    pop_size=1000;
    method='Prop6';
    noise=1;
end
percent7=zeros(size(corrs));
percent8=zeros(size(corrs));
for i=1:size(corrs,2)
    [pt7,pt8,truth]=calcFocals(file,corrs(i),pop_size, method, noise);
    img_idx=@(x) ((imag(x(:,1)))>eps & (imag(x(:,2))))>eps;
    percent7(i)=sum(img_idx(pt7))/pop_size;
    percent8(i)=sum(img_idx(pt8))/pop_size;
end

% histogram errors
figure();
hold on
plot(corrs,percent7);
plot(corrs,percent8);
hold off
xlabel('number of correspondences used')
ylabel('imaginary focal lengths')
legend('7pt','8pt')
saveas(gcf, ['../../results/imaginary_number.jpg']);
saveas(gcf, ['../../results/imaginary_number.fig']);
%toc();
end

function [pt7,pt8,truth]=calcFocals(file,corr,n,method,noise)
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
pt7=zeros(n,2);
pt8=zeros(n,2);
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
    pt7(i,:)=F2f1f2(reshape(Fund,3,3));
    pt7(i,:)=pt7(i,:)*diag([1/A{1}(1) 1/A{2}(1)]).*norm_(i,:);
    
    [Fund,A]=F_features(u1,u2,'Free',testset,3,false);
    pt8(i,:)=F2f1f2(reshape(Fund,3,3));
    pt8(i,:)=pt8(i,:)*diag([1/A{1}(1) 1/A{2}(1)]).*norm_(i,:);
    
    truth(i,:)=truth(i,:).*norm_(i,:);
    rmprintf(repS);
end
end