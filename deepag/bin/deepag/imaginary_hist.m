% A script plots transversal length versus baseline error, also just
% histograms error
% it shows that transversal length has little effect on imaginaryness of
% f-Ratio
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

function [truth,estion,transversal]=imaginary_hist(corr,pop_size, method, noise_out,noise)
if nargin < 1
    %is_ratio=true; %sets whether focal length or ratio is plotted
    noise=1;
    corr=7;
    pop_size=1000;
    method='|F|=0';
    noise_out=0;
end

[truth,estion,transversal]=calcFocals(noise,corr,pop_size, method, noise_out);

% divide data to imaginary and real and calculate f-Ratios
[~,I]=sort(transversal);
transversal=transversal(I);
truth=truth(I,:);
estion=estion(I,:);

real_idx=abs(imag(estion(:,1)))<eps & abs(imag(estion(:,2)))<eps;
imag_idx=abs(imag(estion(:,1)))>eps & abs(imag(estion(:,2)))>eps;
halfreal=not(real_idx | imag_idx);
ratio=abs(estion(:,1)./estion(:,2));
gt_ratio=truth(1,1)./truth(1,2);

if false
% histogram transversal length
figure();
histogram(transversal,20);
xlabel('transversal length');
ylabel('frequency');

% plot relative frequencies of imagninary ratios versus transversal length
figure();
[counts1,centers1]=hist(transversal((real)),20);
[counts2]=hist(transversal(not(real)),centers1);
plot(counts1./counts2);
xlabel('freq_{real}/freq_{imag}');
end 


% plot f-Ratios versus transversal length
figure();
plot(transversal(real_idx),ratio(real_idx),'bx');
hold on
plot(transversal(not(real_idx)),-ratio(not(real_idx)),'rx');
plot([min(transversal) max(transversal)],gt_ratio*ones(1,2));
plot([min(transversal) max(transversal)],-gt_ratio*ones(1,2));
legend({'real f','imaginary f'})
xlabel('transversal length');
ylabel('focal lengths ratio');
axis([0 40 -4 4])
hold off
saveas(gcf, ['../../../results/scatter_imaginary.jpg']);
saveas(gcf, ['../../../results/scatter_imaginary.fig']);


fig=figure();
% histogram errors in ratio for imaginary and real results separately
error_ratio=ratio-gt_ratio;
h1=histogram(error_ratio(real_idx),[-1:0.05:1.5]);
h1.FaceColor='b';
hold on
h2=histogram(error_ratio(imag_idx),[-1:0.05:1.5]); % ,'Normalization','probability'
h2.FaceColor='r';
legend({'real f','imaginary f'})
xlabel('error in ratio');
ylabel('probability');
%title(['error for ' num2str(noise) ' pixel noise']);
saveas(fig, ['../../../results/imaginary_hist_ratio.png']);
saveas(fig, ['../../../results/imaginary_hist_ratio.fig']);

if true
    fig=figure();
    % histogram errors in f for imaginary and real results separately
    error=get_foc_error(estion,truth);
    h1=histogram((error(real_idx)),[0:0.05:1]);
    h1.FaceColor='b';
    hold on
    h2=histogram(error(imag_idx),[0:0.05:1]);
    h2.FaceColor='r';
    legend({'real f','imaginary f'})
    %title(['error for ' num2str(noise) ' pixel noise']);
saveas(gcf, ['../../../results/imaginary_hist.eps'], 'epsc');
saveas(gcf, ['../../../results/imaginary_hist.fig']);
end
end

function [truth,estion,transversal]=calcFocals(noise,corr,n,method,noise_out)
%get n focal length estimations of the bougnoux formula on corr coordinates
%from the generated data

rng(867954151); %the reason I create noise beforehand is that now I can be sure
%for every call of this method noise will be the same
noisy=noise_out*randn(2,n);
noisy(:,1);

estion=zeros(n,2);
support=zeros(n,1);
transversal=zeros(n,1);
truth=repmat([1500 2000],n,1);
for i=1:n
    repS = adprintf({}, [num2str(i), '/', num2str(n)]);
    [F,f1,f2,A, u,P]=getRandomF(noise,60);
    %reshape
    u1=u{1};
    u2=u{2};
    %truncate
    sample=1:corr;
    %sample=randperm(size(u1,2),corr);
    testsample=setdiff(1:size(u1,2),sample);
    testset={u1(:,testsample) u2(:,testsample)};
    u1=u1(:,sample);
    u2=u2(:,sample);
    %noise
    u1(:,1)=u1(:,1)+noisy(:,i);
    
    % find optical axes' transversal length
    [~,~,~,~,~,transversal(i)] = uu2X1([u1;u2],[P{1}; P{2}],'TRAN');
    
    %calculate focal length estimation
    [Fund,A,support(i)]=F_features(u1,u2,method,testset,3,false);
    estion(i,:)=F2f1f2(reshape(Fund,3,3));
    estion(i,:)=estion(i,:)*diag([1/A{1}(1) 1/A{2}(1)]);
    truth(i,:)=truth(i,:);
    rmprintf(repS);
end
end

function [F,f1,f2,A, u,P]=getRandomF(noise,corrnum)
%get random cell array of matrices F 
% F - fundamental matrix
% f1,f2 - focal lengths
% A - scaling matrices
% u - correspondences, beware, has different format

Fparam.f1=1500;
Fparam.f2=2000;
Fparam.per_corr=1;
Fparam.noise=noise;
Fparam.corr=corrnum;

[F,A,u,P]=F_simulate(Fparam,'Free');
f1=Fparam.f1;
f2=Fparam.f2;
end