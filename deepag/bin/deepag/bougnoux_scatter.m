% A script for generating scatter plots of Bougnoux formula results from correspondences file
%
% stat_data:=[absolute_error; a_std; relative_error; r_std]
% data:={error error_proportion} 
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% INRIA, 2016

function [stat_data, data]=bougnoux_scatter(file,corr,pop_size, method)
geometricError=true;

[estion,truth]=calcFocals(file,corr,pop_size, method);
focalscell=num2cell(estion,2);
choice=cellfun(@isreal,focalscell);
unchoice=logical(1-choice);

%calculate statistical data
getprop=@(x) abs(x(:,1))./abs(x(:,2));
divf=@(x,y) abs(x)./abs(y);
diff=@(x,y) abs(x)-abs(y);
if ~geometricError
    divf=diff;
    abs_=@abs;
else
    abs_=@multabs;
end
difg=divf(estion,truth);   
vecprop=@(a,b) deal(divf(a,b),diff(getprop(a),getprop(b)));
[vecA, propA]=vecprop(estion,truth);
[vecR, propR]=vecprop(estion(choice,:),truth(choice,:));
rectify=@(x) reshape(x',[],1);
data{1}=rectify(vecR);
data{2}=rectify(vecA);
data{3}=propR;
data{4}=propA;
%
getmean=@(x)  mean(mean(x,1));
getstd=@(x)  mean(std(x,[],1));
stat_data(:,1)=[abs_(getmean(vecR)) abs_(getmean(vecA))];
stat_data(:,2)=[getstd(vecR) getstd(vecA)];
stat_data(:,3)=[getmean(propR) getmean(propA)];
stat_data(:,4)=[getstd(propR) getstd(propA)];

%plot
scatter(difg(choice,1),difg(choice,2),'b+','DisplayName','real answer');
hold on
if any(unchoice==1)
    scatter(difg(unchoice,1),difg(unchoice,2),'r+','DisplayName','abs(imaginary answer)');
end
%triffles
legend('-DynamicLegend');
title({['Bougnoux formula estimation from ' int2str(corr) ' correspondences.  method = ' method], ...
    [ '[real all] : mean error = ' mat2str(stat_data(:,1),3) '; std = ' mat2str(stat_data(:,2),3)]});
xlabel('abs(f2)');
ylabel('abs(f1)');
method=strrep(method,'|','');
name=['bougnoux/bougnoux_corr=' int2str(corr) '_' method];

%save
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperPosition', [0 0 1024 800]);
saveas(gcf,[name  '.fig']);
saveas(gcf,[name '.jpg']);
hold off
end

function [estion,truth]=calcFocals(file,corr,n,method)
%get n focal length estimations of the bougnoux formula on corr coordinates from the data in
%specified file

%load data
load(file);
u=[corr_tr.u corr_val.u corr_tst.u];
truth=[corr_tr.f corr_val.f corr_tst.f];
norm_=[corr_tr.norm corr_val.norm corr_tst.norm];
%truncate
u=u(:,1:n); truth=truth(:,1:n)'; norm_=norm_(:,1:n)';
estion=zeros(n,2);

for i=1:n
    %reshaping
    uvector=u(:,i);
    points=size(uvector,1)/4;
    u1=reshape(uvector(1:end/2), 2, points);
    u2=reshape(uvector(end/2+1:end), 2, points);
    %truncating
    sample=randperm(size(u1,2),corr);
    u1=u1(:,sample);
    u2=u2(:,sample);
        
    %calculate
    [focal,A]=F_features(u1,u2,method);
    estion(i,:)=F2f1f2(reshape(focal,3,3));
    estion(i,:)=estion(i,:)*diag([1/A{1}(1) 1/A{2}(1)]).*norm_(i,:);
    truth(i,:)=truth(i,:).*norm_(i,:);
end
end

function a=multabs(a)
%absolute value with respect to multiplication. super important function.
if a<1
    a=1/a;
end    
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