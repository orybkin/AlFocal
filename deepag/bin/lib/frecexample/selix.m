% v = selix(m,ix) - select ix in m
% m  = matrix, cell array 
% ix = index 
% v  = m(v)

% T. Pajdla (pajdla@cvut.cz) 2017-04-15
function v = selix(varargin)
v = varargin{1}(varargin{2:end});

