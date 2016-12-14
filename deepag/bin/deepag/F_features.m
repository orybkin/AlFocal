% Calculate feature vector - fundamental matrix - from correspondences
% size(u):=[2 number_of_points]
% (size(u1)==size(u2)):=true
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% INRIA, 2016

function [F,A] = F_features(u1, u2, method, testset,threshold)
if nargin<3
    method='Free';
end
if nargin<4
    testset={NaN};    
end
if nargin<5
    threshold=1/10;
end
%some differents formats of input
if size(u1,1)<3
    u1=[u1; ones(1,size(u1,2))];
    u2=[u2; ones(1,size(u1,2))];
end

if size(testset{1},1)<3 && ~all(any(isnan(testset{1})))
    testset{1}=[testset{1}; ones(1,size(testset{1},2))];
    testset{2}=[testset{2}; ones(1,size(testset{2},2))];
end
%cannot calculate
if size(u1,2)<7 || (size(u1,2)<8 && strcmp(method,'Free'))
    F=NaN;
    A=NaN;
    adprintf({},['they gave me ' num2str(size(u1,2)) ' points']);
    return;
end
% Fundamental matrix
[F,A]  = uu2F({u1,u2},{'None',method});
%[F,A]  = uu2F({u1,u2},{'[-1,1]',method});
if size(F,3)>1
    F=choose_best(F,A,testset,threshold);
end
% normalization
size(F)
size(A{1})
F = inv(A{2})'*F*inv(A{1});
F=reshape(F,9,1);
F = F/norm(F); [~,mi] = max(abs(F)); F = F*sign(F(mi));
end


function F=choose_best(F,A,testset,threshold)
% choice best F via the method of support. Data in 'testset' are used.

anorm=@(x) norm(x(1:2));
if all(any(isnan(testset{1})))
    F=F(:,:,1);
    return;
end
choices=size(F,3);
uset1=testset{1};
uset2=testset{2};
testsize=size(uset1,2);
support=zeros(choices,1);
ll=zeros(testsize,choices);
for i=1:testsize
    u1=uset1(:,i);
    u2=uset2(:,i);
    for j=1:choices
        num=u2'*F(:,:,j)*u1;
        l2=abs(num/anorm(u2'*F(:,:,j)));
        l1=abs(num/anorm(F(:,:,j)*u1));
        l=max(l1,l2);
        ll(i,j)=l
        if l<threshold
            support(j)=support(j)+1;
        end
    end
end
[~,i]=max(support);
F=F(:,:,i);
end