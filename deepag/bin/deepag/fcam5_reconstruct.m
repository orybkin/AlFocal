% Reconstruction 4 focal length computation
%
% Simulation
% Projections
function err=fcam5_reconstruct(K,R,C,X)

P = cellfunuo(@(K,R,C) KRC2P(K,R,C),K,R,C);
u = cellfunuo(@(P) X2u(X,P)+ 0.5*randn(2,size(X,2)),P);
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
% get focal lengths, cameras, reconstruct and select chirality
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
for i=1:numel(u)
    for j=i+1:numel(u)
        % reproject reconstructed points
        ur{i,j}{1} = X2u(Xr{i,j},Pr{i,j}{1});
        ur{i,j}{2} = X2u(Xr{i,j},Pr{i,j}{2});
        % reprojection error
        er{i,j} = [vnorm(u{i}-ur{i,j}{1});vnorm(u{j}-ur{i,j}{2})];
        % find the best transform of the reconstructed points to ground truth points and evaluate the error
        E{i,j} = XY2rts([Xr{i,j};X]);
        XR{i,j} = XE2Y(Xr{i,j},E{i,j}); % points registerd to the ground truth
        EE = [E{i,j}.s*E{i,j}.r E{i,j}.t;0 0 0 1]; % registration transforms
        PR{i,j}{1} = Pr{i,j}{1}*inv(EE); % registered cameras
        PR{i,j}{2} = Pr{i,j}{2}*inv(EE);
        eX{i,j} = XY2rtsErr(E{i,j},[Xr{i,j};X]);
        eXm{i,j} = max(eX{i,j});
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
% Select cameras for the best estimated focal lengths
ix = arrayfun(@(i) argmin(abs(mfe(i)-fe{i})),1:length(mfe)); % best cameras indice for the f's
ci = cell2mat(arrayfun(@(i) fei{i}(:,ix(i)),1:length(ix),'UniformOutput',false)); % best camera pairs 
[ic,jc] = find(ci==ones(size(ci,1),1)*(1:size(ci,2))); % ic is the camera to be is the camera to be used for camera cj in the pair c(:,cj)
Ps = arrayfun(@(i) PR{ci(1,i),ci(2,i)}{ic(i)},1:size(ci,2),'UniformOutput',false); % select cameras
% Initialize points to their geometric median
XX = XR(~cellfun(@isempty,XR)); % all points
nX = size(XX{1},2); % the number of reconstructions 
XX = reshape(cell2mat(XX'),3,nX,[]); % 3D points in a 3D matrix XX(coordinates,points,recontructions)
% get geometric medians
Xs = zeros(3,size(XX,2));
for i=1:size(XX,2)
    x = squeeze(XX(:,i,:)); % all reconstrictions of a point i
    Xs(:,i) = gmedian(x);
end

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
err=eXb;

end
