function [F support] = little_ransac(F,A,testset,threshold,verbose)
% my little_ransac is a small stupid and at the very best partial 
% implementation of this great manifestation of the Monte Carlo 
% methods' genius, RANSAC method.

% choose best F using testset data. A is the scaling, but it is not used.
% what a missed opportunity
global debugg;
if nargin<5
    verbose=false;
end

anorm=@(x) norm(x(1:2));
if all(any(isnan(testset{1})))
    F=F(:,:,1);
    support=NaN;
    return;
end

if size(testset{1},1)<3
    testset{1}=[testset{1}; ones(1,size(testset{1},2))];
    testset{2}=[testset{2}; ones(1,size(testset{2},2))];
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
        ll(i,j)=l;
        if l<threshold
            support(j)=support(j)+1;
        end
    end
end
[~,i]=max(support);
F=F(:,:,i);
if verbose
    support
end
if debugg
    ll;
    size(support);
end
support=support(i)/testsize;
end