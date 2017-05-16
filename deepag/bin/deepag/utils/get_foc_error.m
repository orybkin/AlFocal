% Calculate multiplicative error in focal length
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

function estdata=get_foc_error(estion1,estion2)
% get focal length error
bigger=abs(estion1(:,1))>abs(estion2(:,1));
estdata(bigger)=abs(estion2(bigger,1))./abs(estion1(bigger,1));
estdata(not(bigger))=abs(estion1(not(bigger),1))./abs(estion2(not(bigger),1));
estdata=1-estdata;
end
