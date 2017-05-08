% Calculate feature vector - fundamental matrix - from correspondences
% size(u):=[2 number_of_points]
% (size(u1)==size(u2)):=true
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% INRIA, 2016

function [F,A,support] = F_features(u1, u2, method, testset,threshold,verbose)
if nargin<6
    verbose=false;
end
if nargin<3
    method='Free';
end
if nargin<4
    testset={NaN};    
end
if nargin<5
    threshold=3;
end
%some different formats of input
if size(u1,1)<3
    u1=[u1; ones(1,size(u1,2))];
    u2=[u2; ones(1,size(u1,2))];
end

%cannot calculate
if size(u1,2)<7 || (size(u1,2)<8 && strcmp(method,'Free'))
    F=NaN;
    A=NaN;
    adprintf({},['they gave me ' num2str(size(u1,2)) ' points']);
    return;
end
% Fundamental matrix
[F,A]  = uu2F({u1,u2},{'None',method},testset);
if verbose
    size(F,3)
end
%[F,A]  = uu2F({u1,u2},{'[-1,1]',method});
support=NaN;
if size(F,3)>1
    if all(any(isnan(testset{1}))) % works if testset={nan}. 
        [F]=inliers2F(F,A,{u1 u2});
    else
        [F, support]=voteF(F,A,testset,threshold);
    end
end
% normalization
F = inv(A{2})'*F*inv(A{1});
F=reshape(F,9,1);
F = F/norm(F); [~,mi] = max(abs(F)); F = F*sign(F(mi));
end