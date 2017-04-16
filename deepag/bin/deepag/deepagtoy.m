% deepag.m - Deep Learning 4 Algebraic Geometry - toy example
% T. Pajdla, pajdla@gmail.cz
% O. Rybkin, rybkiole@fel.cvut.cz
% 1 Jan 2016 - 31 Dec 2016
%
%% Initialize
deepagpaths;

%% NEW CODE FROM ZUZANA
% generating
sceneType = {'randombox' 'random'};
pixel = 1/1000;
noise = 1*pixel;  noise=0;
Npoints = 20;
Ncams = 2;
samplesize = 7;
f1 = 1500*pixel;
f2 = 2000*pixel;
fmax=max(f1,f2);
Kgt1 = diag([f1 f1 1]);
Kgt2 = diag([f2 f2 1]);
gtk1 = 0;
gtk2 = 0;
[Pgt M m mgt] = GenerateScene(Npoints, [10*fmax 10*fmax], Ncams, 20*fmax, 30*fmax, 0, noise, [Kgt1;Kgt2], sceneType, [], [], [gtk1,gtk2], true);
Kgt{1} = Kgt1; Kgt{2}=Kgt2;
%ShowCameras(Pgt, Kgt, m, M, true, false, true, 1:7, mgt);
sample=randperm(size(m{1},2),samplesize);
testsample=setdiff(1:Npoints,sample);
u={m{1}(:,sample)/pixel m{2}(:,sample)/pixel};
testset={m{1}(:,testsample)/pixel m{2}(:,testsample)/pixel};
u1=u{1}(:,1);
u{1}(:,1)=u{1}(:,1)+10*randn(2,1);
%% calculating
% 7pt
clear estion3 estion estion6;
[Fsup,Asup]=F_features(u{1},u{2},'|F|=0',testset);
F_t = inv(Asup{2})'*reshape(Fsup,3,3)*inv(Asup{1});
estionsup=F2f1f2(F_t)*diag([1/Asup{1}(1) 1/Asup{2}(1)])


[F3,A]  = uu2F(u,{'None','|F|=0'});
for i=1:size(F3,3)
    F_t = inv(A{2})'*F3(:,:,i)*inv(A{1});
    estion3(i,:)=F2f1f2(F_t)*diag([1/A{1}(1) 1/A{2}(1)]);
end

%% f-Ratio
[Fsup,Asup]=F_features(u{1},u{2},'Prop',testset,5);
F_t = inv(Asup{2})'*reshape(Fsup,3,3)*inv(Asup{1});
estionsup=F2f1f2(F_t)*diag([1/Asup{1}(1) 1/Asup{2}(1)])

[F,A,f]  = uu2F(u,{'None','Prop'},testset);
for i=1:size(F,3)
    F_t = inv(A{2})'*F(:,:,i)*inv(A{1});
    estion(i,:)=F2f1f2(F_t)*diag([1/A{1}(1) 1/A{2}(1)]);
end

%% 6-Ratio
[Fsup,Asup]=F_features(u{1},u{2},'Prop6',testset,3);
F_t = inv(Asup{2})'*reshape(Fsup,3,3)*inv(Asup{1});
estionsup=F2f1f2(F_t)*diag([1/Asup{1}(1) 1/Asup{2}(1)])

[F6,A,f6]  = uu2F(u,{'None','Prop6'},testset,true);
differ=0;
for j=1:size(F6,2)
    for i=1:size(F6{j},3)
        F_t = inv(A{2})'*F6{j}(:,:,i)*inv(A{1});
        estion6{j}(i,:)=F2f1f2(F_t)*diag([1/A{1}(1) 1/A{2}(1)]);
        % differ=differ+norm(estion-estion6{j},'fro');
    end
end
%     for i=1:size(F6,3)
%         F_t = inv(A{2})'*F6(:,:,i)*inv(A{1});
%         estion6(i,:)=F2f1f2(F_t)*diag([1/A{1}(1) 1/A{2}(1)]);
%         % differ=differ+norm(estion-estion6{j},'fro');
%     end

%% plotting
estion3=abs(estion3);
estion=abs(estion);
plot(estion3(:,1),estion3(:,2),'og','DisplayName','7pt solver');
%differ
hold on
plot(estion(:,1),estion(:,2),'or','DisplayName','7pt solver via proportion');
markers={'+' '*' '.' 'x' 's' 'd' 'p' 'h'};
colours={'y' 'm' 'c' 'r' 'g' 'b' 'w' 'k'};
for j=1:size(F6,2)
    estion6{j}=abs(estion6{j});
    [~]=plot(estion6{j}(:,1),estion6{j}(:,2),[markers{j} 'b'],'DisplayName','6pt solver via proportion');
end
hold off
estion3(:,1)./estion3(:,2)
estion(:,1)./estion(:,2)
estion6{1}(:,1)./estion6{1}(:,2)
[~,~]=legend('-DynamicLegend'); % don't change this line - it fixes a Matlab bug



%% OLD CODE ------------------------------------------------------------------------

axdist=100; %pixels
ps.Data    = 'paris'; % Data path
ps.FDemo   = false;    % Cameras & F computation demo
ps.SceneGen= true;   % Scene generator
ps.ptN     = 8;       % The number of points for F computation
ps.ImNS    = 1;     % Image noise ~ N(0,ps.ImNS)
ps.plot = true;
ps.method='|F|=0';
ps.method='Free';
% deepagInit

% focal lengths
f1 = 900;
f2 = 1100;
sc = max(f1,f2);
% 3D points
X=rand(3,15);
X = 5*sc*(X)+repmat([0;0;5*sc],1,size(X,2));
% Internal camera calibration
K1 = [f1  0    0
    0   f1   0
    0   0    1];
K2 = [f2 0  0
    0  f2 0
    0  0  1];
% Rotations
R1 = a2r([1 0 0]',-pi/10);
R2 = a2r([1 0 0]',pi/11);
% Projection matrices
P1 = K1*R1*[eye(3)     [axdist;-sc/2;0]];
P2 = K2*R2*[eye(3)  [0;sc/2;0]];
% Image points
u1 = X2u(X,P1);
u2 = X2u(X,P2);
% add image noise
u1(1:2,:)=u1(1:2,:)+ps.ImNS*randn(size(u1(1:2,:)));
u2(1:2,:)=u2(1:2,:)+ps.ImNS*randn(size(u2(1:2,:)));
% Fundamental matrix
clear e
[F,A]=F_features(u1(:,1:ps.ptN),u2(:,1:ps.ptN),ps.method,{u1(:,ps.ptN:end),u2(:,ps.ptN:end)});
F_t = inv(A{2})'*reshape(F,3,3)*inv(A{1})
Fg = PP2F(P1,P2); % from projection matrices
Fg = Fg/norm(Fg);
festion=F2f1f2(F_t)*diag([1/A{1}(1) 1/A{2}(1)])
%%
if ps.plot
    % Plot cameras & 3D points
    subfig(3,4,1);
    camplot(P1,[],f1/2*[-1 1 1 -1;-1 -1 1 1]);hold; camplot(P2,[],f2/2*[-1 1 1 -1;-1 -1 1 1]);
    plot3d(X,'.');
    axis equal
    title('Cameras & 3D points (''.b'')')
    % Plot images
    subfig(3,4,2);plot3d(u1(1:2,:),'.');axis image;title('Image 1');
    subfig(3,4,3);plot3d(u2(1:2,:),'.');axis image;title('Image 2');
end

if false
    % Plot algebraic error
    subfig(3,4,4);
    e = sum(u2.*(F*u1));
    plot(e,'.-r');axis tight;xlabel('point #');ylabel('error');title('Algebraic error e = u2''*F*u1''');
    % Compare the computed F with the ground truth Fg
    % In image coordinates
    subfig(3,4,5);
    e = [Fg(:)/Fg(3,3)-F(:)/F(3,3)];
    plot(e,'.-r');axis tight;xlabel('F(:)');ylabel('error');title('F(:)-PP2F(P1,P2)(:)');
% Compute the dinfference in normalized image coordinates
Fn = inv(A{2})'*F*inv(A{1});
Fn = Fn/norm(Fn);
Fgn = inv(A{2})'*Fg*inv(A{1});
Fgn = Fgn/norm(Fgn);
% Plot the difference
subfig(3,4,6);
e = [Fgn(:)/norm(Fgn)-Fn(:)/norm(Fn)];
plot(e,'.-r');axis tight;xlabel('Fn(:)');ylabel('error');title('Normalized: F(:)-PP2F(P1,P2)(:)');
% Form the feature vector and the ground truth feature vector
f = Fn(:)/norm(F,'fro'); [~,mi] = max(abs(f)); f = f*sign(f(mi));
fg = Fgn(:)/norm(Fg,'fro'); [~,mi] = max(abs(fg)); fg = fg*sign(fg(mi));
subfig(3,4,7);
plot(1:length(f),f+1,'.-',1:length(fg),fg+1,'.-g');axis tight;xlabel('element');ylabel('value');title('b - f, g - fg');
% Bougnoux formula aplied on normalized data and recomputed to original data
ff = F2f1f2(F)
fn = F2f1f2(Fn);
fgn = F2f1f2(Fgn);
fnb = fn*diag([1/A{1}(1) 1/A{2}(1)])
fgnb = fgn*diag([1/A{1}(1) 1/A{2}(1)]);

end
if ps.plot && false
    % plot estimated f's
    plot([1 2],[[f1 f2];ff;fnb]','.','markersize',20);
    xlabel('focal length index'); ylabel('pixels');
    title(sprintf('b-truth, g-est, r-norm im(f) = %.1f %.1f',imag(fnb(1)),imag(fnb(2))));
end
