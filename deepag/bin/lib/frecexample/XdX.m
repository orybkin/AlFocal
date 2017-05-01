% d = XdX(X) - kD point L2 distance matrix
% 
% X ... k x n point coordinates
% d ... n x n point symmetric matrix of L2 distances

% (c) T.Pajdla, pajdla@cvut.cz, 2017-04-14
function d = XdX(X)
if nargin>0
    % 2D
    % dm = ([repmat(x,[size(x,2),1])-repmat(x(:),[1,size(x,2)])]);
    % d  = sqrt(dm(1:2:end,:).^2+dm(2:2:end,:).^2);
    % kD
    d = reshape(vnorm([kron(X,ones(1,size(X,2)))-kron(ones(1,size(X,2)),X)]),[size(X,2) size(X,2)]);
else % unit tests 
    X = [1 3 5;2 4 6];
    d0 = [0      2.8284 5.6569
          2.8284 0      2.8284
          5.6569 2.8284 0     ];
    d = all(all(abs(d0-XdX(X))<1e-4));
end