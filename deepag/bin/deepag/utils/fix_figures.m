% mass creation of pdf from fig. useful e.g. when switching to another matlab
% version
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% INRIA, 2016
stdir=cd;
%cd bougnoux;
cd /media/intelligent-agent/Sky/Cloud/Opera/CMP/FocalNet/FocalNet/results;
figs=dir('*.fig');
for i=1:size(figs,1)
    fig=openfig(figs(i).name,'invisible');
    saveas(fig,[figs(i).name(1:end-3) 'eps'],'epsc');
end