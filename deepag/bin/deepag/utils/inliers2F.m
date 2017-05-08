function [F] = voteF(F,A,testset)
% choose best F with sampson error. A is the scaling, but it is not used.
% what a missed opportunity
global debugg;
anorm=@(x) norm(x(1:2));
if iscell(F)
    F=cell2mat(reshape(F,1,1,[]));
end


if size(testset{1},1)<3
    testset{1}=[testset{1}; ones(1,size(testset{1},2))];
    testset{2}=[testset{2}; ones(1,size(testset{2},2))];
end

choices=size(F,3);
uset1=testset{1};
uset2=testset{2};
for i=1:choices
    err(i)=sampson_error(F(:,:,i),uset1,uset2);
end
[~,i]=min(err);
F=F(:,:,i);
end