% [i,m] = argmin(A,[],dim,nanflag) - argument of min
% A,[],dim,nanflag = input as for min
% [i,m] = [m,i] = min(A,[],dim,nanflag)

% T. Pajdla (pajdla@cvut.cz) 2017-04-15
function [i,m] = argmin(varargin)
if nargin>0
    [m,i] = min(varargin{:});
else
    i = true;
end

