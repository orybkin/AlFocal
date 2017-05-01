% Reconstruction 4 focal length computation
%
close all
clear all
% Simulation
% Scene
X = [0 1 1 0 0 1 2 0 1 2
     0 0 1 1 0 0 1 1 1 2
     0 0 0 0 1 2 1 2 1 2];
X = X + 0.1*randn(size(X));
X = [X 2*rand(3,40)+[0;0;0]*ones(1,40)];
% Cameras
% 1
f0{1} = 1000;
K{1} = [f0{1} 0 500;0 f0{1} 500;0 0 1];
R{1} = a2r([1;10;0],20/180*pi);
C{1} = 3*[2;1;-2];
% 2
f0{2} = 2000;
K{2} = [f0{2} 0 500;0 f0{2} 500;0 0 1];     
R{2} = a2r([7;10;2],-20/180*pi);
C{2} = 3*[-1;1;-2];
% 3
f0{3} = 3000;
K{3} = [f0{3} 0 500;0 f0{3} 500;0 0 1];     
R{3} = a2r([10;5;1],30/180*pi);
C{3} = 3*[1;-1;-2];
if true  
    % 4
    f0{4} = 4000;
    K{4} = [f0{4} 0 500;0 f0{4} 500;0 0 1];
    R{4} = a2r([10;-10;1],-10/180*pi);
    C{4} = 5*[0.5;0.5;-2];
    % 5
    f0{5} = 2000;
    K{5} = [f0{5} 0 500;0 f0{5} 500;0 0 1];
    R{5} = a2r([-3;45;-3],-75/180*pi);
    C{5} = 3*[-2;0;0];
end
P = cellfunuo(@(K,R,C) KRC2P(K,R,C),K,R,C);
f1 = subfig(3,4,1); plot3d(X,'.k');hold; cellfunuo(@(P) camplot(P,1),P); grid; axis equal; view([0 -80]);
title('Scene & cameras'); xlabel('x'); ylabel('y'); zlabel('z');
% Projections
u = cellfunuo(@(P) X2u(X,P)+ 0.5*randn(2,size(X,2)),P);
%
clr = {'b' 'g' 'r' 'k' 'm' 'c' 'y'};
f2 = subfig(3,4,2); hold;
cellfunuo(@(u,c) plot3d(u,['.' c]),u,clr(1:numel(u))); 
axis equal;
axis([0 1000 0 1000]);
% Fundamental matrix
for i=1:numel(u)
    for j=i+1:numel(u)
        F{i,j} = uu2F({a2h(u{i}),a2h(u{j})},{'[-1,1]','|F|=0'});
        for k=1:size(F{i,j},3)
            F{i,j}(:,:,k) = F{i,j}(:,:,k)/norm(F{i,j}(:,:,k));
            e{i,j}(k) = max(EGeomErr(F{i,j}(:,:,k),[a2h(u{i});a2h(u{j})]));
        end
        % select the right F
        [e{i,j},ix] = min(e{i,j});
        F{i,j} = F{i,j}(:,:,ix);
    end
end
title(['Images & epipolar lines, max(e) = ' sprintf('%.2g ',[e{:}])]);
% epipolar lines
for i=1:numel(u)
    for j=i+1:numel(u)
        li= F{i,j}'*a2h(u{j});
        lj = F{i,j}*a2h(u{i});
        set(plotline(li),'color',clr{j});
        set(plotline(lj),'color',clr{i});
    end
end
% get focal lengths, cameras, reconstruct and select
PR = {}; Pr = {}; XR = {}; Xr = {}; v = {}; z = {}; e = {}; a = {};
for i=1:numel(u)
    for j=i+1:numel(u)
        % principal points
        p{i} = K{i}(1:2,3);
        p{j} = K{j}(1:2,3);
        % Recover f1, f2 from Ef
        f{i,j} = F2f1f2(F{i,j},{p{i},p{j}},'Bougnoux');
        % Internal calibration
        Kr1 = [f{i,j}(1) 0 p{i}(1);0 f{i,j}(1) p{i}(2);0 0 1];
        Kr2 = [f{i,j}(2) 0 p{j}(1);0 f{i,j}(2) p{j}(2);0 0 1];
        % Essential matrix
        E{i,j} = Kr2'*F{i,j}*Kr1;
        % Reconstructed calibrated projection matrices
        Pc{i,j} = E2PP(E{i,j});
        % Reconstruct and select the best reconstruction = most points in front of cameras        
        for k=1:numel(Pc{i,j}{2})
            % add internal calibration to Pr's
            Pr{i,j}{1}{k} = Kr1*Pc{i,j}{1}{k};
            Pr{i,j}{2}{k} = Kr2*Pc{i,j}{2}{k};
            % reconstruct
            [Xr{i,j}{k},v{i,j}{k},z{i,j}{k},e{i,j}{k},a{i,j}{k}] = uu2X([u{i};u{j}],[Pr{i,j}{1}{k};Pr{i,j}{2}{k}]);
        end
        % number of points in front of both cameras
        np = cell2mat(cellfunuo(@(z) sum(all(z>0)),z{i,j}));
        [np,ix] = max(np);
        Xr{i,j} = Xr{i,j}{ix};
        Pr{i,j}{1} = Pr{i,j}{1}{ix};
        Pr{i,j}{2} = Pr{i,j}{2}{ix};        
    end
end
% evaluate the initial recosntruction
f3 = subfig(3,4,3); hold;
for i=1:numel(u)
    for j=i+1:numel(u)
        % reproject reconstructed points
        ur{i,j}{1} = X2u(Xr{i,j},Pr{i,j}{1});
        ur{i,j}{2} = X2u(Xr{i,j},Pr{i,j}{2});
        % reprojection error
        er{i,j} = [vnorm(u{i}-ur{i,j}{1});vnorm(u{j}-ur{i,j}{2})];
        figure(f2); plot3d(ur{i,j}{1},['o' clr{i}]);plot3d(ur{i,j}{2},['o' clr{j}]);
        % find the best transform of the reconstructed points to ground truth points and evaluate the error
        E{i,j} = XY2rts([Xr{i,j};X]);
        XR{i,j} = XE2Y(Xr{i,j},E{i,j}); % points registerd to the ground truth
        EE = [E{i,j}.s*E{i,j}.r E{i,j}.t;0 0 0 1]; % registration transforms
        PR{i,j}{1} = Pr{i,j}{1}*inv(EE); % registered cameras
        PR{i,j}{2} = Pr{i,j}{2}*inv(EE);
        eX{i,j} = XY2rtsErr(E{i,j},[Xr{i,j};X]);
        eXm{i,j} = max(eX{i,j});
        figure(f3); plot(eX{i,j},'color',clr{i});
    end
end
axis tight;
title(['Point registration error, max=' sprintf('%.2g ',[eXm{:}])]);xlabel('point');ylabel('error');
% show cameras
figure(f1);
for i=1:numel(u)
    for j=i+1:numel(u)
        plot3d(XR{i,j},'.b');
        set(camplot(PR{i,j}{1},1),'color',clr{i});set(camplot(PR{i,j}{2},1),'color',clr{j});
    end
end
% Extract focal length estimates
fe  = cell(numel(u),0);
fei = cell(numel(u),0);
for i=1:numel(u)
    fe{i} = [];
    fei{i} = [];
end
for i=1:numel(u)
    for j=i+1:numel(u)
        fe{i} = [fe{i} f{i,j}(1)]; % focal length
        fei{i} = [fei{i} [i;j]]; % the pair of cameras leading to the focal length
        fe{j} = [fe{j} f{i,j}(2)]; % focal length
        fei{j} = [fei{j} [i;j]]; % the pair of cameras leading to the focal length
    end
end
% Estimate the focal length as the median of all estimated f's
mfe = cell2mat(cellfunuo(@(x) median(x),fe));
% Show focal langth estimates
f4 = subfig(3,4,4); hold
arrayfun(@(i) plot(repmat(i,1,numel(fe{i})),fe{i},['.' clr{i}]),1:numel(fe));
arrayfun(@(i) plot(i,f0{i},['s' clr{i}]),1:numel(fe));
arrayfun(@(i) plot(i,mfe(i),['o' clr{i}]),1:numel(fe));
axis('tight');xlabel('camera');ylabel('f');
% Select cameras for the best estimated focal lengths
ix = arrayfun(@(i) argmin(abs(mfe(i)-fe{i})),1:length(mfe)); % best cameras indice for the f's
ci = cell2mat(arrayfun(@(i) fei{i}(:,ix(i)),1:length(ix),'UniformOutput',false)); % best camera pairs 
[ic,jc] = find(ci==ones(size(ci,1),1)*(1:size(ci,2))); % ic is the camera to be is the camera to be used for camera cj in the pair c(:,cj)
Ps = arrayfun(@(i) PR{ci(1,i),ci(2,i)}{ic(i)},1:size(ci,2),'UniformOutput',false); % select cameras
% Initialize points to their geoetric median
XX = XR(~cellfun(@isempty,XR)); % all points
nX = size(XX{1},2); % the number of reconstructions 
XX = reshape(cell2mat(XX'),3,nX,[]); % 3D points in a 3D matrix XX(coordinates,points,recontructions)
% get geometric medians
Xs = zeros(3,size(XX,2));
for i=1:size(XX,2)
    x = squeeze(XX(:,i,:)); % all reconstrictions of a point i
    Xs(:,i) = gmedian(x);
end
f5 = subfig(2,3,4);
plot3d(X,'.k'); hold; plot3d(Xs,'.b'); 
cellfunuo(@(P) camplot(P,1),P); 
arrayfun(@(i) set(camplot(Ps{i},1),'color',clr{i}),1:numel(Ps));
grid; axis equal; view([0 -80]);
 
% Bundle adjustment
opBA.proj_func = 1; % projection function ar,s,x0,y0 is fixed, other adjusted
opBA.robustify = 0; % robust loss function (0 - least squares, 1 - Huber's loss), default 0
opBA.maxerr = 3; % maximal residual after BA for inliers
opBA.tolerance = 1e-12; % stopping criterion - change of error, parameters and ...
opBA.p_tolerance = 1e-12; % stopping criterion - change of optimized variables
opBA.g_tolerance = 1e-12; % stopping criterion - size of the gradient
opBA.f_tolerance = 1e-12; % stopping criterion - change of the objective
opBA.num_iterations = 100; % the maximum number of LM iterations, default 5
opBA.num_threads = 1; % number of threads, default 1
opBA.verbose = 0; % info printout on
opBA.constant_cameras = 0; % first camera fixed
opBA.constant_points = 0; % no fixed cameras
% do BA
[ub,Pb,Xb,eb] = uPXBA(u,Ps,Xs,opBA);
% error on points
Eb = XY2rts([Xb;X]);
eXb = XY2rtsErr(Eb,[Xb;X]);
figure(f3); plot(eXb,'.-g','linewidth',3);
% show cameras
EE = [Eb.s*Eb.r Eb.t;0 0 0 1];
figure(f5); plot3d(h2a(EE*a2h(Xb)),'og');
arrayfun(@(i) set(camplot(Pb{i},1),'color',clr{i},'linewidth',2),1:numel(Pb));
title(sprintf('Rec & -BA & +BA, eX = %.2g',mean(eXb))); xlabel('x'); ylabel('y'); zlabel('z');
% show final focal lengths
fb = cell2mat(cellfunuo(@(P) mean(selix(P2KRC(P),[1 5])),Pb));
figure(f4)
arrayfun(@(i) plot(i,fb(i),['*' clr{i}]),1:length(fb));
title('s GT, . Pairs, o -BA, * +BA');