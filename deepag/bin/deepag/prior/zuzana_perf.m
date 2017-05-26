% The script histograms Prior solver results versus baseline
%
% Oleh Rybkin, rybkiole@fel.cvut.cz
% CMP, 2017

function [algs,prior]=noise_test(corr,pop_size, cnoise, fnoises)
if nargin < 1
    fnoises=[0 1 5 10 20 30 50]./100;
    cnoise = [0];
    corr=7.;
    pop_size=50;
end

addpath('./utils')
addpath('./hartley')
addpath('./solver')
addpath('./scene')
addpath('./common')

f1gt = 1700;
f2gt = 2500;
u1gt = 20;
u2gt = 20;
v1gt = 10;
v2gt = 10;

pixel = 1/1000; 

[algs,prior]=calcFocals(cnoise,fnoises,corr,pop_size, f1gt, f2gt, u1gt, v1gt, u2gt, v2gt);

% histogram errors
%figure();
%hold on
%histhist(sort(zuzana.iter),20,1,'-g');
%histhist(sort(hartley.iter),20,1,'-r');
%hold off
%legend('zuzana iterations', 'hartley iterations')

%----vykreslenie

cfgstyles = {'-' '-.*' '--' ':'};
methodsNames = { 'New' 'Hartley' 'Bougnoux' 'Hartley ratio'};
name = '';   
cfgs = [1];
methods = [1 2 4 3]
cols = {'g' 'b' 'r'  'k' 'c' 'm' 'y'};       % colors

f1err = true;

if f1err
    
    collector = [];
    ncnt = length(fnoises);
    
    %  allmes{cfg}{focal}{alg}(noise, measurement)
    i = 1;
    ii = 1;
    
    
    hnd = figure;
    axes_handle=axes('fontsize', 16);
    hold on;
    set(hnd, 'Name', name);
    title(name);
    lims =  [0   1.3];
    ylim(axes_handle, lims);
    xlim(axes_handle, [0.4 length(fnoises)+0.5]);
    
    xt = [];
    xn = {};
    for ns=1:length(fnoises)
        xt(ns) = ns*1;
        xn{ns} = num2str(fnoises(ns));
        plot(-1,-1, cols{ns});
    end
    
    set(gca,'XTick', xt);
    set(gca,'XTickLabel', xn);
    
    xlabel('noise in prior', 'fontsize', 16);
    ylabel('Relative error of estimated f_1', 'fontsize', 16);
    
    for ns=2:length(fnoises)
        plot([ns-0.5 ns-0.5], [lims(1) lims(2)], ['-'], 'LineWidth', 1, 'color', [0.9 0.9 0.9]);
    end
    
    metshift = 0.9/length(methods);
    xsize = metshift/2-0.02;
    xxsize = metshift/4-0.02;
    ofs = (1 + xsize) - (metshift * length(methods)) / 2;
    
    cfg = cfgs(1);
    
    %plot([0 10], [1 1], 'Color', 'cyan');
    
    scollector = [];
    for j=methods
        
        xstart = ofs;
        xpos = xstart;
        
        for ns=1:length(fnoises)
            
            x = [];
                      
                mx =abs(algs{j}(ns).f(:,1)-(f1gt*pixel))/(f1gt*pixel);
                x = [x; mx(:)];
                
            
            
            [outlier,loadj,upadj,yy,q1,q3,n2,med,n1] = mb_boxutil(x, 1,1.5,0);
            
            plot([xpos xpos], [q1 n2], ['--' cols{ii}], 'LineWidth', 1);
            plot([xpos xpos], [q3 n1], ['--' cols{ii}], 'LineWidth', 1);
            
            plot(xpos * ones(1, length(yy)), yy, ['+' cols{ii}], 'MarkerSize', 3);
            
            plot([xpos-xxsize xpos+xxsize], [q1 q1], ['-' cols{ii}], 'LineWidth', 1);
            plot([xpos-xxsize xpos+xxsize], [q3 q3], ['-' cols{ii}], 'LineWidth', 1);
            
            plot([xpos-xsize xpos+xsize xpos+xsize xpos-xsize xpos-xsize], [n1 n1 n2 n2 n1], ['-' cols{ii}], 'LineWidth', 1);
            
            plot([xpos-xsize xpos+xsize], [med med], ['-' cols{ii}], 'LineWidth', 2);
            
            stat = [n2, med, n1];
            
            %scollector(:, i) = stat;
            %collector(:, i) = x;
            i=i+1;
            xpos = xpos + 1;
        end
        
        %X = (1:size(lambdas,2)) + xstart -1;
  
        
        ofs = ofs + metshift;
        ii = ii + 1;
        
        legend(methodsNames(methods));
        
    end
end


f2err = true;

if f2err
    
    collector = [];
    ncnt = length(fnoises);
    
    %  allmes{cfg}{focal}{alg}(noise, measurement)
    i = 1;
    ii = 1;
    
    
    hnd = figure;
    axes_handle=axes('fontsize', 16);
    hold on;
    set(hnd, 'Name', name);
    title(name);
    lims =  [0   1.3];
    ylim(axes_handle, lims);
    xlim(axes_handle, [0.4 length(fnoises)+0.5]);
    
    xt = [];
    xn = {};
    for ns=1:length(fnoises)
        xt(ns) = ns*1;
        xn{ns} = num2str(fnoises(ns));
        plot(-1,-1, cols{ns});
    end
    
    set(gca,'XTick', xt);
    set(gca,'XTickLabel', xn);
    
    xlabel('noise in prior', 'fontsize', 16);
    ylabel('Relative error of estimated f_2', 'fontsize', 16);
    
    for ns=2:length(fnoises)
        plot([ns-0.5 ns-0.5], [lims(1) lims(2)], ['-'], 'LineWidth', 1, 'color', [0.9 0.9 0.9]);
    end
    
    metshift = 0.9/length(methods);
    xsize = metshift/2-0.02;
    xxsize = metshift/4-0.02;
    ofs = (1 + xsize) - (metshift * length(methods)) / 2;
    
    cfg = cfgs(1);
    
    %plot([0 10], [1 1], 'Color', 'cyan');
    
    scollector = [];
    for j=methods
        
        xstart = ofs;
        xpos = xstart;
        
        for ns=1:length(fnoises)
            
            x = [];
                      
                mx =abs(algs{j}(ns).f(:,2)-(f2gt*pixel))/(f2gt*pixel);
                x = [x; mx(:)];
                
            
            
            [outlier,loadj,upadj,yy,q1,q3,n2,med,n1] = mb_boxutil(x, 1,1.5,0);
            
            plot([xpos xpos], [q1 n2], ['--' cols{ii}], 'LineWidth', 1);
            plot([xpos xpos], [q3 n1], ['--' cols{ii}], 'LineWidth', 1);
            
            plot(xpos * ones(1, length(yy)), yy, ['+' cols{ii}], 'MarkerSize', 3);
            
            plot([xpos-xxsize xpos+xxsize], [q1 q1], ['-' cols{ii}], 'LineWidth', 1);
            plot([xpos-xxsize xpos+xxsize], [q3 q3], ['-' cols{ii}], 'LineWidth', 1);
            
            plot([xpos-xsize xpos+xsize xpos+xsize xpos-xsize xpos-xsize], [n1 n1 n2 n2 n1], ['-' cols{ii}], 'LineWidth', 1);
            
            plot([xpos-xsize xpos+xsize], [med med], ['-' cols{ii}], 'LineWidth', 2);
            
            stat = [n2, med, n1];
            
            %scollector(:, i) = stat;
            %collector(:, i) = x;
            i=i+1;
            xpos = xpos + 1;
        end
        
        %X = (1:size(lambdas,2)) + xstart -1;
  
        
        ofs = ofs + metshift;
        ii = ii + 1;
        
        legend(methodsNames(methods));
        
    end
end

methods = [1,2, 4];
p1err = true;

if p1err
    
    collector = [];
    ncnt = length(fnoises);
    
    %  allmes{cfg}{focal}{alg}(noise, measurement)
    i = 1;
    ii = 1;
    
    
    hnd = figure;
    axes_handle=axes('fontsize', 16);
    hold on;
    set(hnd, 'Name', name);
    title(name);
    lims =  [0   0.5];
    ylim(axes_handle, lims);
    xlim(axes_handle, [0.4 length(fnoises)+0.5]);
    
    xt = [];
    xn = {};
    for ns=1:length(fnoises)
        xt(ns) = ns*1;
        xn{ns} = num2str(fnoises(ns));
        plot(-1,-1, cols{ns});
    end
    
    set(gca,'XTick', xt);
    set(gca,'XTickLabel', xn);
    
    xlabel('noise in prior', 'fontsize', 16);
    ylabel('Error of principal point 1', 'fontsize', 16);
    
    for ns=2:length(fnoises)
        plot([ns-0.5 ns-0.5], [lims(1) lims(2)], ['-'], 'LineWidth', 1, 'color', [0.9 0.9 0.9]);
    end
    
    metshift = 0.9/length(methods);
    xsize = metshift/2-0.02;
    xxsize = metshift/4-0.02;
    ofs = (1 + xsize) - (metshift * length(methods)) / 2;
    
    cfg = cfgs(1);
    
    %plot([0 10], [1 1], 'Color', 'cyan');
    
    scollector = [];
    for j=methods
        
        xstart = ofs;
        xpos = xstart;
        
        for ns=1:length(fnoises)
            
            x = [];
                      
                mx =sqrt(sum(((algs{j}(ns).p1-([u1gt,v1gt]*pixel)).^2)')')
                x = [x; mx(:)];
                
            
            
            [outlier,loadj,upadj,yy,q1,q3,n2,med,n1] = mb_boxutil(x, 1,1.5,0);
            
            plot([xpos xpos], [q1 n2], ['--' cols{ii}], 'LineWidth', 1);
            plot([xpos xpos], [q3 n1], ['--' cols{ii}], 'LineWidth', 1);
            
            plot(xpos * ones(1, length(yy)), yy, ['+' cols{ii}], 'MarkerSize', 3);
            
            plot([xpos-xxsize xpos+xxsize], [q1 q1], ['-' cols{ii}], 'LineWidth', 1);
            plot([xpos-xxsize xpos+xxsize], [q3 q3], ['-' cols{ii}], 'LineWidth', 1);
            
            plot([xpos-xsize xpos+xsize xpos+xsize xpos-xsize xpos-xsize], [n1 n1 n2 n2 n1], ['-' cols{ii}], 'LineWidth', 1);
            
            plot([xpos-xsize xpos+xsize], [med med], ['-' cols{ii}], 'LineWidth', 2);
            
            stat = [n2, med, n1];
            
            %scollector(:, i) = stat;
            %collector(:, i) = x;
            i=i+1;
            xpos = xpos + 1;
        end
        
        %X = (1:size(lambdas,2)) + xstart -1;
  
        
        ofs = ofs + metshift;
        ii = ii + 1;
        
        legend(methodsNames(methods));
        
    end
end

methods = [1,2, 4];
p2err = true;

if p2err
    
    collector = [];
    ncnt = length(fnoises);
    
    %  allmes{cfg}{focal}{alg}(noise, measurement)
    i = 1;
    ii = 1;
    
    
    hnd = figure;
    axes_handle=axes('fontsize', 16);
    hold on;
    set(hnd, 'Name', name);
    title(name);
    lims =  [0   0.5];
    ylim(axes_handle, lims);
    xlim(axes_handle, [0.4 length(fnoises)+0.5]);
    
    xt = [];
    xn = {};
    for ns=1:length(fnoises)
        xt(ns) = ns*1;
        xn{ns} = num2str(fnoises(ns));
        plot(-1,-1, cols{ns});
    end
    
    set(gca,'XTick', xt);
    set(gca,'XTickLabel', xn);
    
    xlabel('noise in prior', 'fontsize', 16);
    ylabel('Error of principal point 2', 'fontsize', 16);
    
    for ns=2:length(fnoises)
        plot([ns-0.5 ns-0.5], [lims(1) lims(2)], ['-'], 'LineWidth', 1, 'color', [0.9 0.9 0.9]);
    end
    
    metshift = 0.9/length(methods);
    xsize = metshift/2-0.02;
    xxsize = metshift/4-0.02;
    ofs = (1 + xsize) - (metshift * length(methods)) / 2;
    
    cfg = cfgs(1);
    
    %plot([0 10], [1 1], 'Color', 'cyan');
    
    scollector = [];
    for j=methods
        
        xstart = ofs;
        xpos = xstart;
        
        for ns=1:length(fnoises)
            
            x = [];
                      
                mx =sqrt(sum(((algs{j}(ns).p2-([u2gt,v2gt]*pixel)).^2)')')
                x = [x; mx(:)];
                
            
            
            [outlier,loadj,upadj,yy,q1,q3,n2,med,n1] = mb_boxutil(x, 1,1.5,0);
            
            plot([xpos xpos], [q1 n2], ['--' cols{ii}], 'LineWidth', 1);
            plot([xpos xpos], [q3 n1], ['--' cols{ii}], 'LineWidth', 1);
            
            plot(xpos * ones(1, length(yy)), yy, ['+' cols{ii}], 'MarkerSize', 3);
            
            plot([xpos-xxsize xpos+xxsize], [q1 q1], ['-' cols{ii}], 'LineWidth', 1);
            plot([xpos-xxsize xpos+xxsize], [q3 q3], ['-' cols{ii}], 'LineWidth', 1);
            
            plot([xpos-xsize xpos+xsize xpos+xsize xpos-xsize xpos-xsize], [n1 n1 n2 n2 n1], ['-' cols{ii}], 'LineWidth', 1);
            
            plot([xpos-xsize xpos+xsize], [med med], ['-' cols{ii}], 'LineWidth', 2);
            
            stat = [n2, med, n1];
            
            %scollector(:, i) = stat;
            %collector(:, i) = x;
            i=i+1;
            xpos = xpos + 1;
        end
        
        %X = (1:size(lambdas,2)) + xstart -1;
  
        
        ofs = ofs + metshift;
        ii = ii + 1;
        
        legend(methodsNames(methods));
        
    end
end


methods = [1 2 4]
niter = 1;
if niter
    
    collector = [];
    ncnt = length(fnoises);
    
    %  allmes{cfg}{focal}{alg}(noise, measurement)
    i = 1;
    ii = 1;
    
    
    hnd = figure;
    axes_handle=axes('fontsize', 16);
    hold on;
    set(hnd, 'Name', name);
    title(name);
    lims =  [0   190];
    ylim(axes_handle, lims);
    xlim(axes_handle, [0.4 length(fnoises)+0.5]);
    
    xt = [];
    xn = {};
    for ns=1:length(fnoises)
        xt(ns) = ns*1;
        xn{ns} = num2str(fnoises(ns));
        plot(-1,-1, cols{ns});
    end
    
    set(gca,'XTick', xt);
    set(gca,'XTickLabel', xn);
    
    xlabel('noise in prior', 'fontsize', 16);
    ylabel('Iterations', 'fontsize', 16);
    
    for ns=2:length(fnoises)
        plot([ns-0.5 ns-0.5], [lims(1) lims(2)], ['-'], 'LineWidth', 1, 'color', [0.9 0.9 0.9]);
    end
    
    metshift = 0.9/length(methods);
    xsize = metshift/2-0.02;
    xxsize = metshift/4-0.02;
    ofs = (1 + xsize) - (metshift * length(methods)) / 2;
    
    cfg = cfgs(1);
    
    %plot([0 10], [1 1], 'Color', 'cyan');
    
    scollector = [];
    for j=methods
        
        xstart = ofs;
        xpos = xstart;
        
        for ns=1:length(fnoises)
            
            x = [];
                      
                mx =algs{j}(ns).iter(:);
                x = [x; mx(:)];
                
            
            
            [outlier,loadj,upadj,yy,q1,q3,n2,med,n1] = mb_boxutil(x, 1,1.5,0);
            
            plot([xpos xpos], [q1 n2], ['--' cols{ii}], 'LineWidth', 1);
            plot([xpos xpos], [q3 n1], ['--' cols{ii}], 'LineWidth', 1);
            
            plot(xpos * ones(1, length(yy)), yy, ['+' cols{ii}], 'MarkerSize', 3);
            
            plot([xpos-xxsize xpos+xxsize], [q1 q1], ['-' cols{ii}], 'LineWidth', 1);
            plot([xpos-xxsize xpos+xxsize], [q3 q3], ['-' cols{ii}], 'LineWidth', 1);
            
            plot([xpos-xsize xpos+xsize xpos+xsize xpos-xsize xpos-xsize], [n1 n1 n2 n2 n1], ['-' cols{ii}], 'LineWidth', 1);
            
            plot([xpos-xsize xpos+xsize], [med med], ['-' cols{ii}], 'LineWidth', 2);
            
            stat = [n2, med, n1];
            
            %scollector(:, i) = stat;
            %collector(:, i) = x;
            i=i+1;
            xpos = xpos + 1;
        end
        
        %X = (1:size(lambdas,2)) + xstart -1;
  
        
        ofs = ofs + metshift;
        ii = ii + 1;
        
        legend(methodsNames(methods));
        
    end
end
    
end

%toc();


function [algs,prior]=calcFocals(cnoise, fnoises, corr,n, f1gt, f2gt, u1gt, v1gt, u2gt, v2gt)
global debugg;
%get n focal length estimations of the bougnoux formula on corr coordinates from the data in
%specified file



%rng(867954152); 
rng(156350212); 
rng(34135153);

for j = 1:size(fnoises,2)
    
    
    
prior(j).f=zeros(n,2);
prior(j).p1=zeros(n,2);
prior(j).p1=zeros(n,2);

algs{1}(j).f=zeros(n,2);
algs{1}(j).p1=zeros(n,2);
algs{1}(j).p2=zeros(n,2);
algs{1}(j).iter=zeros(n,1);
algs{1}(j).converge=zeros(n,1);

algs{2}(j).f=zeros(n,2);
algs{2}(j).p1=zeros(n,2);
algs{2}(j).p2=zeros(n,2);
algs{2}(j).iter=zeros(n,1);
algs{2}(j).converge=zeros(n,1);

algs{3}(j).f=zeros(n,2);

for i=1:n
    repS = adprintf({}, [num2str(i), '/', num2str(n)]);
    %generate
    sceneType = {'randombox' 'random'};
    pixel = 1/1000;
    noise = cnoise*pixel; 
    Npoints = 20;
    Ncams = 2;
    f1 = f1gt*pixel;
    f2 = f2gt*pixel;
    p1= [u1gt*pixel v1gt*pixel];
    p2= [u2gt*pixel v2gt*pixel];
    fmax=max(f1,f2);
    Kgt1 = diag([f1 f1 1]);
    Kgt2 = diag([f2 f2 1]);
    Kgt1(1:2,3)=p1;
    Kgt2(1:2,3)=p2;
    gtk1 = 0;
    gtk2 = 0;
    [Pgt M m mgt] = GenerateScene(Npoints, [10*fmax 10*fmax], Ncams, 10*fmax, 30*fmax, 0, noise, [Kgt1;Kgt2], sceneType, [], [], [gtk1,gtk2], true);
    Kgt{1} = Kgt1; Kgt{2}=Kgt2;
    %ShowCameras(Pgt, Kgt, m, M, true, false, true, 1:7, mgt);
    sample=randperm(size(m{1},2),corr);
    testsample=setdiff(1:Npoints,sample);
    u1=m{1}(:,sample);
    u2=m{2}(:,sample);
    testset={m{1}(:,testsample) m{2}(:,testsample)};
    
    
    %calculate baseline estimation
    [F,A]=F_features(u1,u2,'|F|=0',testset,0.001);
    F=reshape(F,3,3);

    algs{3}(j).f(i,:)=F2f1f2(F);
    algs{3}(j).f(i,:)=algs{3}(j).f(i,:)*diag([1/A{1}(1) 1/A{2}(1)]);
    
        
    ratio=abs(F2ratio(F));
    real_idx=@(x) ((imag(x(:,1))<eps) & (imag(x(:,2)))<eps);
    img_idx=@(x) ((imag(x(:,1))>eps) & (imag(x(:,2)))>eps);
    if ~real_idx(algs{3}(j).f(i,:)) && ~img_idx(algs{3}(j).f(i,:))
        ratio=1;
    end
    
    
    %f1=baseline(i,1);
    %f2=baseline(i,2);
    %p1=[0;0];
    %p2=[0;0];
    
    % Zuzana
    %f1prior =f1+ 1.5*rand(1);
    %u1prior = p1(1)+0.1*rand(1);
    %v1prior = p1(2)+0.1*rand(1);
    
    f1prior =f1+ fnoises(j)*f1;
    u1prior = p1(1)+fnoises(j)*p1(1);
    v1prior = p1(2)+fnoises(j)*p1(2);
    
    u1prior = 0;
    v1prior = 0;
    
    w1 = 0.1;
    w2 = 1;
    
   % f2prior = f2+ 1.5*rand(1);
   % u2prior = p2(1)+0.05*rand(1);
   % v2prior = p2(2)+0.05*rand(1);
   
    f2prior =f2+ fnoises(j)*f2;
    u2prior = p2(1)+fnoises(j)*p2(1);
    v2prior = p2(2)+fnoises(j)*p2(2);
    
    u2prior = 0;
    v2prior = 0;
    
    
    prior(j).f(i,:)=[f1prior f2prior];
    prior(j).p1(i,:)=[u1prior v1prior];
    prior(j).p2(i,:)=[u2prior v2prior];
    w3 =0.1;
    w4 =1;
   
    tic()
    [f1_zuz, u1_zuz, v1_zuz, f2_zuz, u2_zuz, v2_zuz, l1, l2, err, iter] = f1_f2_from_F(F, f1prior, u1prior, v1prior,f2prior, u2prior, v2prior, w1, w2, w3, w4);
    toc()
    
    algs{1}(j).converge(i)=iter>150;
    algs{1}(j).iter(i)=iter;
    algs{1}(j).f(i,:)=[f1_zuz f2_zuz];
    algs{1}(j).p1(i,:)=[u1_zuz v1_zuz];
    algs{1}(j).p2(i,:)=[u2_zuz v2_zuz];
    
    % Hartley
     tic()
     [F, p_hartley, f_hartley, output] = uu2F_hartley(F,u1,u2,[f1prior; f2prior],{[u1prior; v1prior] [u1prior; v1prior]});
     toc()
     algs{2}(j).converge(i)=output.exitflag;
     algs{2}(j).iter(i)=output.iterations;
     algs{2}(j).f(i,:)=f_hartley;
     algs{2}(j).p1(i,:)=p_hartley{1}';
     algs{2}(j).p2(i,:)=p_hartley{2}';
     %norm(f_hartley-[f1 f2])
          
     
     % Hartley prior
     tic()
     [F2, p_hartley2, f_hartley2, output2] = uu2F_hartley(F,u1,u2,[f1prior; f2prior],{[u1prior; v1prior] [u1prior; v1prior]},0,ratio);
     toc()
     algs{4}(j).converge(i)=output2.exitflag;
     algs{4}(j).iter(i)=output2.iterations;
     algs{4}(j).f(i,:)=f_hartley2;
     algs{4}(j).p1(i,:)=p_hartley2{1}';
     algs{4}(j).p2(i,:)=p_hartley2{2}';
     %norm(f_hartley-[f1 f2])
    rmprintf(repS);
end
end
end