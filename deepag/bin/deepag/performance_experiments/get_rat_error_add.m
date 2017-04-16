% Calculate additive error in f-ratio
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

function estdata=get_rat_error_add(estion,truth)
% get additive focal lengths ratio error
estion=estion(:,1)./estion(:,2);
truth=truth(:,1)./truth(:,2);
estdata=abs(estion)-abs(truth);
% remove outliers
estdata=estdata(estdata<10);
end