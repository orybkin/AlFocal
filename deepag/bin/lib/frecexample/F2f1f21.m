% f = F2f1f2(F,[p,mth]) - Focal lengths from Fundamental matrix
%
% F   = Fundamental matrix with diagonal K1 = [f1 0 p1(1); 0 f2 p1(2); 0 0 1], 
%                                        K2 = [f2 0 p1(2); 0 f2 p2(2); 0 0 1]
% p   = {p1 p2} princial points ({[0;0] [0;0]} implicit)
% mth = method
%      'Bougnoux' - Bougnoux formula from HZ (implicit), doesnot work for some pure translations
%      'Robust-3' - Switches between 3 formulas to cover pure translations  
%
% T. Pajdla, pajdla@cvut.cz, 2016-09-30
function f = F2f1f2(F,p,mth)
if nargin>0
    if nargin<3
        mth = 'Bougnoux';
    end
    if nargin<2 || isempty(p)
        p = {[0;0] [0;0]};
    end
    switch mth
        case 'Bougnoux'
            % principal points
            p1 = a2h(p{1});
            p2 = a2h(p{2});
            % epipoles
            [e2,~,e1] = svd(F);
            e1 = e1(:,3);
            e2 = e2(:,3);
            % the formula
            II = diag([1 1 0]);
            f(1) = sqrt(abs((-p2'*xx(e2)*II*F *(p1*p1')*F'*p2)/(p2'*xx(e2)*II*F *II*F'*p2)));
            f(2) = sqrt(abs((-p1'*xx(e1)*II*F'*(p2*p2')*F *p1)/(p1'*xx(e1)*II*F'*II*F *p1)));
        case 'Bougnoux-N'
            A1 = [1 0 -p{1}(1);0 1 -p{1}(2);0 0 1]; 
            A2 = [1 0 -p{2}(1);0 1 -p{2}(2);0 0 1]; 
            F = inv(A2)'*F*inv(A1);
            % principal points in [0,0]
            p1 = A1*a2h(p{1});
            p2 = A2*a2h(p{2});
            % epipoles
            [e2,~,e1] = svd(F);
            e1 = e1(:,3);
            e2 = e2(:,3);
            % the formula
            II = diag([1 1 0]);
            f(1) = sqrt(abs((-p2'*xx(e2)*II*F *(p1*p1')*F'*p2)/(p2'*xx(e2)*II*F *II*F'*p2)));
            f(2) = sqrt(abs((-p1'*xx(e1)*II*F'*(p2*p2')*F *p1)/(p1'*xx(e1)*II*F'*II*F *p1)));   
        case 'Robust-3'
            % move the principal point to [0,0]
            A1 = [1 0 -p{1}(1);0 1 -p{1}(2);0 0 1]; 
            A2 = [1 0 -p{2}(1);0 1 -p{2}(2);0 0 1]; 
            F = inv(A2)'*F*inv(A1);
            s = [min(svd(F(:,[1 2]))) min(svd(F(:,[2 3]))) min(svd(F(:,[3 1])))]; % numerical ranks of two-column submatrices
            [sm,ix] = max(s); % the most non-singular submatrix
            if sm<eps, error('F has rank 1'); end
            ix = 0;
            switch ix
                case 0
                    f(1) = sqrt((-F(1,1)*F(2,3)*F(3,1)*F(3,3) - F(1,2)*F(2,3)*F(3,2)*F(3,3) + F(1,3)*F(2,1)*F(3,1)*F(3,3) + F(1,3)*F(2,2)*F(3,2)*F(3,3))/...
                                (F(1,1)^2*F(1,3)*F(2,3) - F(1,1)*F(1,3)^2*F(2,1) + F(1,1)*F(2,1)*F(2,3)^2 + F(1,2)^2*F(1,3)*F(2,3) - F(1,2)*F(1,3)^2*F(2,2) + F(1,2)*F(2,2)*F(2,3)^2 - F(1,3)*F(2,1)^2*F(2,3) - F(1,3)*F(2,2)^2*F(2,3))); 
                    f(2) = sqrt((-F(1,1)*F(3,2)*F(1,3)*F(3,3) - F(1,2)*F(3,2)*F(2,3)*F(3,3) + F(3,1)*F(1,2)*F(3,1)*F(3,3) + F(3,1)*F(2,2)*F(2,3)*F(3,3))/...
                                (F(1,1)^2*F(3,1)*F(3,2) - F(1,1)*F(3,1)^2*F(1,2) + F(1,1)*F(1,2)*F(3,2)^2 + F(2,1)^2*F(3,1)*F(3,2) - F(2,1)*F(3,1)^2*F(2,2) + F(2,1)*F(2,2)*F(3,2)^2 - F(3,1)*F(1,2)^2*F(3,2) - F(3,1)*F(2,2)^2*F(3,2)));   
                case 1
                     % nm = - f11*f23*f31*f33 - f12*f23*f32*f33 + f13*f21*f31*f33 + f13*f22*f32*f33
                     %         2                 2                 2      2                 2                 2          2             2
                     % dn = f11 f13*f23 - f11*f13 f21 + f11*f21*f23  + f12 f13*f23 - f12*f13 f22 + f12*f22*f23  - f13*f21 f23 - f13*f22 f23 
                     %      - f11*f13*f32*f33 + f12*f13*f31*f33 - f21*f23*f32*f33 + f22*f23*f31*f33
                     % nm = - f11*f13*f32*f33 - f12*f13*f31*f33 + f21*f23*f32*f33 - f22*f23*f31*f33
                     %         2                     2              2      2             2                     2              2      2
                     % dn = f11 f31*f32 - f11*f12*f31  + f11*f12*f32  - f12 f31*f32 + f21 f31*f32 - f21*f22*f31  + f21*f22*f32  - f22 f31*f32
                     f(1) = sqrt((F(1,2)*F(1,3)*F(3,3)^2-F(1,3)^2*F(3,2)*F(3,3)+F(2,2)*F(2,3)*F(3,3)^2-F(2,3)^2*F(3,2)*F(3,3)/...
                                 (F(1,1)*F(1,2)*F(3,1)*F(3,3)-F(1,1)*F(1,3)*F(3,1)*F(3,2)+F(1,2)^2*F(3,2)*F(3,3)-F(1,2)*F(1,3)*F(3,2)^2+F(2,1)*F(2,2)*F(3,1)*F(3,3)-F(2,1)*F(2,3)*F(3,1)*F(3,2)+F(2,2)^2*F(3,2)*F(3,3)-F(2,2)*F(2,3)*F(3,2)^2)));
                     f(2) = sqrt((F(2,1)*F(3,1)*F(3,3)^2+F(2,2)*F(3,2)*F(3,3)^2-F(2,3)*F(3,1)^2*F(3,3)-F(2,3)*F(3,2)^2*F(3,3))/...
                                 (F(1,1)*F(1,3)*F(2,1)*F(3,3)-F(1,1)*F(1,3)*F(2,3)*F(3,1)*+F(1,2)*F(1,3)*F(2,2)*F(3,3)-F(1,2)*F(1,3)*F(2,3)*F(3,2)+F(2,1)^2*F(2,3)*F(3,3)-F(2,1)*F(2,3)^2*F(3,1)+F(2,2)^2*F(2,3)*F(3,3)-F(2,2)*F(2,3)^2*F(3,2)));
                case 2
                     f(1) = sqrt((F(1,1)*F(1,3)*F(3,3)^2-F(1,3)^2*F(3,1)*F(3,3)+F(2,1)*F(2,3)*F(3,3)^2-F(2,3)^2*F(3,1)*F(3,3)/...
                                 (F(1,1)^2*F(3,1)*F(3,3)+F(1,1)*F(1,2)*F(3,2)*F(3,3)-F(1,1)*F(1,3)*F(3,1)^2-F(1,2)*F(1,3)*F(3,1)*F(3,2)+F(2,1)^2*F(3,1)*F(3,3)+F(2,1)*F(2,2)*F(3,2)*F(3,3)-F(2,1)*F(2,3)*F(3,1)^2-F(2,2)*F(2,3)*F(3,1)*F(3,2))));
                     f(2) = sqrt((F(1,1)*F(3,1)*F(3,3)^2+F(1,2)*F(3,2)*F(3,3)^2-F(1,3)*F(3,1)^2*F(3,3)-F(1,3)*F(3,2)^2*F(3,3))/...
                                 (F(1,1)^2*F(1,3)*F(3,3)-F(1,1)*F(1,3)^2*F(3,1)+F(1,1)*F(2,1)*F(2,3)*F(3,3)+F(1,2)^2*F(1,3)*F(3,3)-F(1,2)*F(1,3)^2*F(3,2)+F(1,2)*F(2,2)*F(2,3)*F(3,3)-F(1,3)*F(2,1)*F(2,3)*F(3,1)-F(1,3)*F(2,2)*F(2,3)*F(3,2)));
                case 3
                     f(1) = sqrt((F(1,1)*F(1,3)*F(3,2)*F(3,3)-F(1,2)*F(1,3)*F(3,1)*F(3,3)+F(2,1)*F(2,3)*F(3,2)*F(3,3)-F(2,2)*F(2,3)*F(3,1)*F(3,3))/...
                                 (F(1,1)^2*F(3,1)*F(3,2)-F(1,1)*F(1,2)*F(3,1)^2+F(1,1)*F(1,2)*F(3,2)^2-F(1,2)^2*F(3,1)*F(3,2)+F(2,1)^2*F(3,1)*F(3,2)-F(2,1)*F(2,2)*F(3,1)^2+F(2,1)*F(2,2)*F(3,2)^2-F(2,2)^2*F(3,1)*F(3,2)));
                     f(2) = sqrt((F(1,1)*F(3,1)*F(2,3)*F(3,3)-F(2,1)*F(3,1)*F(1,3)*F(3,3)+F(1,2)*F(3,2)*F(2,3)*F(3,3)-F(2,2)*F(3,2)*F(1,3)*F(3,3))/...
                                 (F(1,1)^2*F(1,3)*F(2,3)-F(1,1)*F(2,1)*F(1,3)^2+F(1,1)*F(2,1)*F(2,3)^2-F(2,1)^2*F(1,3)*F(2,3)+F(1,2)^2*F(1,3)*F(2,3)-F(1,2)*F(2,2)*F(1,3)^2+F(1,2)*F(2,2)*F(2,3)^2-F(2,2)^2*F(1,3)*F(2,3)));
            end
        otherwise
            f = [NaN NaN];
    end
else % unit tests
    p = {[500;600] [500;400]};
    fg = [1100 900];
    K1 = [fg(1)  0     p{1}(1)
          0      fg(1) p{1}(2)
          0      0     1];
    K2 = [fg(2) 0     p{2}(1)
          0     fg(2) p{2}(2)
          0     0     1];
    % Rotations
    R1 = eye(3);
    R2 = a2r(rand(3,1),pi/10);
    % Projection matrices
    P1 = K1*R1*[eye(3)       [-1000;0;0]];
    P2 = K2*R2*[eye(3)  1000*[rand(2,1);0]];
    F = PP2F(P1,P2);
    % Test 1: Bougnoux formula    
    ff = F2f1f2(F,p,'Bougnoux');
    f(1) = all(abs(ff-fg)<1e-8);
    % Test 2: Bougnoux formula + normalization   
    ff = F2f1f2(F,p,'Bougnoux-N');
    f(2) = all(abs(ff-fg)<1e-8);
    % Test 3: Robust-3
    ff = F2f1f2(F,p,'Robust-3');
    f(3) = all(abs(ff-fg)<1e-8);
end
