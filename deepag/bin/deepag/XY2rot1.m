% R = XY2rot(X1,X2) - optimal rotation registration
%
% R = arg min ||X2 - r*X1||
%
% X1 ... 3 x n points
% X2 ... 3 x n points
% R  ... 3 x 3 rotation

% T. Pajdla, pajdla@cmp.felk.cvut.cz
% 2015-08-09
function R = XY2rot1(X1,X2)
l2R = @(l) [cos(l) -sin(l) ; sin(l) cos(l)];
if nargin>0
    n=size(X1,1);
    % direction correlation matrix
    K = zeros(n,n);				
    for i=1:size(X1,2)						
      K = K + X2(:,i)*X1(:,i)';	
    end			
    [U,~,V] = svd(K);
    R = U*V';
    if n==2
        l=atan(sum(X1(1,:).*X2(2,:)-X1(2,:).*X2(1,:))/sum(X1(1,:).*X2(1,:)+X1(2,:).*X2(2,:)));
        R=l2R(l-pi);
    end
else
    X1=rand(2,7);
    X1=X1-mean(X1);
    R=l2R(2)
    X2=R*X1;
    XY2rot1(X1,X2)
end