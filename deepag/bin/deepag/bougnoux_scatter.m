% A script for generating scatter plots of Bougnoux formula results from correspondences file
%
% stat_data:=[absolute_error; a_std; relative_error; r_std]
% data:={absolute_error a_std relative_error r_std
%        absolute_error_img a_std_img relative_error_img r_std_img}
%
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% INRIA, 2016

function [stat_data, data]=bougnoux_scatter(file,corr,pop_size, method, noise,plotting)
if strcmp(method,'Free') && corr<8
    data={[] [] [] []}';
    stat_data=nan(2,4);
    return;
end
geometricError=true;
[estion,truth]=calcFocals(file,corr,pop_size, method, noise);
estion(:,1)./estion(:,2);
%truth(:,1)./truth(:,2)

tic();
focalscell=num2cell(estion,2);
real_idx=cellfun(@isreal,focalscell);
img_idx=logical(1-real_idx);

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
prop_errs=divf(estion,truth);   
%vec - vector, the focal lengths themselves, prop - their proportion f1/f2
vec_prop=@(a,b) deal(divf(a,b),diff(getprop(a),getprop(b)));
[vecImg, propImg]=vec_prop(estion(img_idx,:),truth(img_idx,:));
[vecReal, propReal]=vec_prop(estion(real_idx,:),truth(real_idx,:));
rectify=@(x) reshape(x',[],1);
data{1}=rectify(vecReal);
data{2}=rectify(vecImg);
data{3}=propReal;
data{4}=propImg;
%
getmean=@(x)  mean(mean(x,1));
getstd=@(x)  mean(std(x,[],1));
stat_data(:,1)=[abs_(getmean(vecReal)) abs_(getmean(vecImg))];
stat_data(:,2)=[getstd(vecReal) getstd(vecImg)];
stat_data(:,3)=[getmean(propReal) getmean(propImg)];
stat_data(:,4)=[getstd(propReal) getstd(propImg)];
%toc();
%plot
if plotting
    scatter(abs(estion(real_idx,1)),abs(estion(real_idx,2)),'b+','DisplayName','real focal');
    hold on
    if any(img_idx==1)
        scatter(abs(estion(img_idx,1)),abs(estion(img_idx,2)),'r+','DisplayName','abs(imaginary focal)');
    end
    scatter(1500,2000,40,'go','DisplayName','ground truth');
    plot([0 3000],[0 4000],'y','DisplayName','correct proportion line');
    axis([0 4000 0 4000]);
    %triffles
    [~,~]=legend('-DynamicLegend'); % don't change this line - it fixes a Matlab bug
    title({['Bougnoux formula estimation from ' int2str(corr) ' correspondences.  method = ' method], ...
        [ '[real all] : mean error = ' mat2str(stat_data(:,1),3) '; std = ' mat2str(stat_data(:,2),3)]});
    xlabel('abs(f2)');
    ylabel('abs(f1)');
    method=strrep(method,'|','');
    name=['bougnoux/scatter/bougnoux_corr=' int2str(corr) '_' method];
    saveas(gcf,[name  '.fig']);
    saveas(gcf,[name '.jpg']);
    hold off
end
%toc();
end

function [estion,truth]=calcFocals(file,corr,n,method,noise)
global debugg;
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
rng(867954152); %the reason I create noise beforehand is that now I can be sure 
%for every call of this method the noise will be the same, which allows for
%comparison
noisy=noise*randn(2,n);
noisy(:,1);


for i=1:n
    %reshape
    uvector=u(:,i);
    points=size(uvector,1)/4;
    u1=reshape(uvector(1:end/2), 2, points);
    u2=reshape(uvector(end/2+1:end), 2, points);
    %truncate
    sample=1:corr;
    %sample=randperm(size(u1,2),corr);
    testsample=setdiff(1:size(u1,2),sample);
    testset={u1(:,testsample) u2(:,testsample)};
    u1=u1(:,sample);
    u2=u2(:,sample);
    %noise
    u1(:,1)=u1(:,1)+noisy(:,i);
        
 
    crucial=5;
    %calculate
    if i==crucial 
        debugg=0;
        
        method;
    else
        debugg=0;
    end
        
    [Fund,A]=F_features(u1,u2,method,testset,3,false);
    %method
    %little_ransac(reshape(Fund,3,3),A,testset,3,false);
    estion(i,:)=F2f1f2(reshape(Fund,3,3));
    estion(i,:)=estion(i,:)*diag([1/A{1}(1) 1/A{2}(1)]).*norm_(i,:);
    truth(i,:)=truth(i,:).*norm_(i,:);
    if i==crucial
       % estion(i,:)
     %   estion(i,1)/estion(i,2)
    end
    
end
end

function a=multabs(a)
%absolute value with respect to multiplication. super important function.
if abs(a)<1
    a=1/a;
end    
end

function Fparam=getFparam(corr)
% DATED
%get a structure with camera system parameters
Fparam.corr=corr; % number of correspondences
Fparam.f1=900;
Fparam.f2=1100;
Fparam.points=rand(3,Fparam.corr);
Fparam.raxis=rand(3,1); % r
Fparam.tpos=rand(2,1);
end