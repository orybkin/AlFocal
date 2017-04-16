% The script plots focal length errors against optical axes distance 
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

function [estion,truth]=degeneracy_perf(noise,corr,pop_size, method)
global parallel_axes;
if nargin < 1
    parallel_axes=true;
    corr=7;
    pop_size=100;
    method='|F|=0';
    noise=[1 1e-1 1e-2 1e-3];
    %noise=[1e-1];
    if parallel_axes
        % elevation angle
        param=[0 0.1 0.2 0.5 1 2 5 10];
    else
        % optical axes distance
        param=[0 0.1 0.2 0.5 1 2 5 10 50];
    end
    %axdist=[0 0.2 1 5];
    colors={'-' '.-' '-' '.-' '-' '.-' '-' '.-' '-' '.-' '-' '.-'};
end

params=size(param,2);
noises=size(noise,2);
N=floor(sqrt(noises));
M=ceil(noises/N);
[estion,truth]=calcFocals(corr,pop_size, method, noise,param);

% plot focal errors
figure();
for j=1:noises
    subplot(N,M,j);
    for i=1:params
        cumhist(sort(get_foc_error(estion(:,:,i,j),truth)),20,2,[colors{i}]);
    end
    title(['noise ' num2str(noise(j)) 'px']);
    xlabel('error magnitude');
    ylabel('frequency, %');
end
if parallel_axes
    legend(arrayfun(@(x) {[num2str(x) char(176)]},param))
else
    legend(arrayfun(@(x) {[num2str(x) ' cm']},param))
end
%suptitle('errors in focal length w.r.t. axes distance')

% plot ratio errors
figure();
for j=1:noises
    subplot(N,M,j);
    for i=1:params
        cumhist(sort(get_rat_error_mult(estion(:,:,i,j),truth)),20,2,[colors{i}]);
    end
    title(['noise ' num2str(noise(j)) 'px']);
    xlabel('error magnitude');
    ylabel('frequency, %');
end
legend(arrayfun(@(x) {num2str(x)},param))
suptitle('multiplicative errors in ratio w.r.t. axes distance')

% plot ratio errors
figure();
for j=1:noises
    subplot(N,M,j);
    for i=1:params
        cumhist(sort(get_rat_error_add(estion(:,:,i,j),truth)),20,2,[colors{i}]);
    end
    title(['noise ' num2str(noise(j)) 'px']);
    xlabel('error magnitude');
    ylabel('frequency, %');
    xlim([-1,1]);
end
legend(arrayfun(@(x) {num2str(x)},param))
suptitle('additive errors in ratio w.r.t. axes distance')

end

function [estion,truth]=calcFocals(corr,n,method,noise,axdist)
%get n focal length estimations of the bougnoux formula on corr coordinates from the data in
%specified file

axdists=size(axdist,2);
noises=size(noise,2);
estion=zeros(n,2,axdists,noises);
truth=zeros(n,2);

rng(867954152);

for i=1:n
    repS = adprintf({}, [num2str(i), '/', num2str(n)]);
    % SIMULATE DATA
    
    % focal lengths
    truth(i,:)=[900, 1100];
    [X, K]=getPoints(truth);
    for j=1:axdists
        for k=1:noises
            [~,estion(i,:,j,k)]=getF(X,K,noise(k),method,corr,axdist(j));
        end
    end
    
    % what to do with imaginary things
    rmprintf(repS);
end
end

function [X, K] = getPoints(f)
% 3D points
X=rand(3,15);
X = 500*(X)+repmat([0;0;500],1,size(X,2));
% Internal camera calibration
K1 = [f(1)  0    0
    0   f(1)   0
    0   0    1];
K2 = [f(2) 0  0
    0  f(2) 0
    0  0  1];
K={K1,K2};
end

function [F, f]=getF(X,K,noise,method,points,axdist)
global parallel_axes;
% the distances are measured in cm
if parallel_axes
    % Rotations
    R1 = a2r([0 1 0]',degtorad(axdist));
    R2 = a2r([0 1 0]',0);
    R={R1,R2};
    % translations
    t1=[0;-100/2;0];
    t2=[0;100/2;0];
else
    % Rotations
    R1 = a2r([1 0 0]',-pi/10);
    R2 = a2r([1 0 0]',pi/11);
    R={R1,R2};
    % translations
    t1=[axdist;-100/2;0];
    t2=[0;100/2;0];
end
% Projection matrices
P1 = K{1}*R{1}*[eye(3)  t1];
P2 = K{2}*R{2}*[eye(3)  t2];
% Image points
u1 = X2u(X,P1);
u2 = X2u(X,P2);
% add image noise
u1(1:2,:)=u1(1:2,:)+noise*randn(size(u1(1:2,:)));
u2(1:2,:)=u2(1:2,:)+noise*randn(size(u2(1:2,:)));
% Fundamental matrix
[F,A]=F_features(u1(:,1:points),u2(:,1:points),method,{u1(:,points:end),u2(:,points:end)});
F_t = inv(A{2})'*reshape(F,3,3)*inv(A{1});
f=F2f1f2(F_t)*diag([1/A{1}(1) 1/A{2}(1)]);
end