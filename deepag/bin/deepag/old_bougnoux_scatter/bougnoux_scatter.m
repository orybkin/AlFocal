% A script for generating scatter plots of Bougnoux formula results
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% INRIA, 2016

function stat_data=bougnoux_scatter(Fparam, runs, method)
%stat_data:=[absolute_error; a_std; relative_error; r_std]
two_focals=false;

focals=calcFocals(runs,Fparam, method);

focalscell=num2cell(focals,2);
choice=cellfun(@isreal,focalscell);
unchoice=logical(1-choice);

%calculate statistical data
a=focals(choice,:); b=focals;
if two_focals
    getaerr=@(x)  (mean(abs(x),1)-[Fparam.f1 Fparam.f2]);
    getastd=@(x)  (std(abs(x),[],1));
    stat_data(:,1)=[getaerr(a)];
    stat_data(:,2)=[getastd(a) ];
    getrerr=@(x)  mean(abs(x(:,1))./abs(x(:,2)),1)-Fparam.f1/Fparam.f2;
    getrstd=@(x)  std(abs(x(:,1))./abs(x(:,2)),[],1);
    stat_data(:,3)=[getrerr(a) getrerr(b)];
    stat_data(:,4)=[getrstd(a) getrstd(b)];
else
    getaerr=@(x)  mean(mean(abs(x),1)-[Fparam.f1 Fparam.f2]);
    getastd=@(x)  mean(std(abs(x),[],1));
    stat_data(:,1)=[getaerr(a) getaerr(b)];
    stat_data(:,2)=[getastd(a) getastd(b)];
    getrerr=@(x)  mean(abs(x(:,1))./abs(x(:,2)),1)-Fparam.f1/Fparam.f2;
    getrstd=@(x)  std(abs(x(:,1))./abs(x(:,2)),[],1);
    stat_data(:,3)=[getrerr(a) getrerr(b)];
    stat_data(:,4)=[getrstd(a) getrstd(b)];
end

%plot
scatter(abs(focals(choice,1)),abs(focals(choice,2)),'b+','DisplayName','real focal');
hold on
if any(unchoice==1)
    scatter(abs(focals(unchoice,1)),abs(focals(unchoice,2)),'r+','DisplayName','abs(imaginary focal)');
end
scatter(900,1100,40,'go','DisplayName','ground truth');
plot([0 3600],[0 4400],'y','DisplayName','correct proportion line');

%
legend('-DynamicLegend');
title({['Bougnoux formula estimation. ' mat2str(Fparam.corr) ' correspondences, ' mat2str(runs) ' runs. method = ' method], ...
    [ '[real all] : mean error = ' mat2str(stat_data(:,1),3) '; std = ' mat2str(stat_data(:,2),3)]});
xlabel('abs(f2)');
ylabel('abs(f1)');
axis([0 3600 0 4400]);
name=['bougnoux/bougnoux_' int2str(Fparam.corr) 'cor_' int2str(runs) 'runs_' method];

%save
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [0 0 1024 800]);
saveas(gcf,[name  '.fig']);
saveas(gcf,[name '.jpg']);
hold off
end

function fgnb=calcFocals(n,Fparam,method)
%get n focal length estimations from the bougnoux formula
fgnb=zeros(n,2);
tic();
for i=1:n
    Fparam.points=rand(3,Fparam.corr);
    temp=FDemo(Fparam,1,method);
    fgnb(i,:)=temp;
end
toc();
end

function Fparam=getFparam(corr)
%get a structure with camera system parameters
Fparam.corr=corr; % number of correspondences
Fparam.f1=900;
Fparam.f2=1100;
Fparam.points=rand(3,Fparam.corr);
Fparam.raxis=rand(3,1); % r
Fparam.tpos=rand(2,1);
end
