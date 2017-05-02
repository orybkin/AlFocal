% Calculate multiplicative error in focal length
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

function estdata=get_foc_error(estion,truth)
% get focal length error
bigger=abs(estion(:,1))>abs(truth(:,1));
estdata(bigger)=abs(truth(bigger,1))./abs(estion(bigger,1));
estdata(not(bigger))=abs(estion(not(bigger),1))./abs(truth(not(bigger),1));
estdata=1-estdata;
end
