function [F support] = voteF(F,A,testset,threshold,verbose)
% support method for chosing F

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
unique(support);
if size(unique(support),1)==1
    ll;
    warning(['Support couldnt decide. See "voteF". Support is ' mat2str(support) '. Distances are ' mat2str(ll)]);
end
[~,i]=max(support);
F=F(:,:,i);
if verbose
    support;
end
if debugg
    ll;
    size(support);
end
support=support(i)/testsize;
end