% Calculate multiplicative error in f-ratio
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

function estdata=get_rat_error_mult(estion,truth)
estion=estion(:,1)./estion(:,2);
truth=truth(:,1)./truth(:,2);

estdata=get_foc_error(estion,truth);
end
