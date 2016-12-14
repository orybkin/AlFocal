function [eq, known, unknown, kngroups, cfg, algB] = sw6pt()
%%
% 6point focal length problem
% http://cmp.felk.cvut.cz/minimal/6_pt_relative.php

%%
setpaths;

%%
% create symbolic matrices 
% these matrices represent basis of null space of linearized equation y'*E*x = 0
% matrices F1, F2, F3 will be treated as known
F1 = gbs_Matrix('F1%d%d', 3 ,3);
F2 = gbs_Matrix('F2%d%d', 3 ,3);
F3 = gbs_Matrix('F3%d%d', 3 ,3);

% create fundamental matrix as a linear combination of the null space,
% w = 1/focal^2
% prop = f2/f1
% Q = f^-2 * K^2
% (symbolicaly - variables x, y, w are unknown)
syms x y w prop real;
F = x*F1 + y*F2 + F3;

% build equations
%  1. fundamental matrix is singular
eq(1) = det(F); 

Q1 = diag([1 1 w]);
Q2 = diag([1 1 w*prop^2]);

% calibrate fundamental matrix using unknown focal length and use "trace"
% constraint for obtained essential matrix
Ft = transpose(F);
te = 2*(F*Q2*Ft*Q1)*F - trace(F*Q2*Ft*Q1)*F;
eq(2:10) = te(:);

%%
% collect known & unknown variables
unknown = {'x' 'y' 'w'};
vars = transpose([F1(:); F2(:); F3(:)]);
known = {'prop'};
for var = vars
    known = [known {char(var)}];
end

%%

% define variable groups (optional)
kngroups = ones(9,1)*[1 2 3];
kngroups = [4; kngroups(:)];
kngroups=[];

% define configuration
cfg = gbs_InitConfig();

% no algB yet computed
algB = [];

end