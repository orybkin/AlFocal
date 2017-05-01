% f = F2f1f2(F,[p,mth]) - Focal lengths from Fundamental matrix
%
% F   = Fundamental matrix with diagonal K1 = diag([f1 f2 1]), K2 = diag([f2 f2 1])
% p   = {p1 p2} princial points ({[0;0] [0;0]} implicit)
% mth = method
%      'Bougnoux' - Bougnoux formula from HZ (implicit)
%
% T. Pajdla, pajdla@cvut.cz, 2016-08-28
function f_ratio = F2ratio(F,p,mth)
if nargin>0
    if nargin<3
        mth = 'Polynomial';
    end
    if nargin<2
        p = {[0;0] [0;0]};
    end
    switch mth
        case 'Bougnoux'
            % epipoles
            F=F/F(3,3);
            if false
                [e2,~,e1] = svd(F);
                e1 = e1(:,3);
                e2 = e2(:,3);
            else
                e1=cross(F(1,:),F(2,:))';
                e2=cross(F(:,1),F(:,2));
            end
            e1 = e1/e1(3);
            e2 = e2/e2(3);
            
            II = diag([1 1 0]);
            p1 = a2h(p{1});
            p2 = a2h(p{2});
            f_ratio(1) = sqrt((-p2'*xx(e2)*II*F *(p1*p1')*F'*p2)/(p2'*xx(e2)*II*F *II*F'*p2));
            f_ratio(2) = sqrt((-p1'*xx(e1)*II*F'*(p2*p2')*F *p1)/(p1'*xx(e1)*II*F'*II*F *p1));
            if any(~isreal(f_ratio)) && false
                [(-p1'*xx(e1)*II*F'*(p2*p2')*F *p1) (p1'*xx(e1)*II*F'*II*F *p1)]
                [(-p2'*xx(e2)*II*F *(p1*p1')*F'*p2) (p2'*xx(e2)*II*F *II*F'*p2)]
                1;
            end
            f_ratio=f_ratio(2)/f_ratio(1);
            % formulae from Macaulay - see documentation (unfortunately it does not exist)
        case 'Polynomial1'
            f_ratio(1)=-F(3,3)*(F(1,2)*F(1,3)*F(3,3)-F(1,3)^2*F(3,2)+F(2,2)*F(2,3)*F(3,3)-F(2,3)^2*F(3,2))/ ...
                (F(1,1)*F(1,2)*F(3,1)*F(3,3)-F(1,1)*F(1,3)*F(3,1)*F(3,2)+F(1,2)^2*F(3,2)*F(3,3)-F(1,2)*F(1,3)*F(3,2)^2+F(2,1)*F(2,2)*F(3,1)*F(3,3)-F(2,1)*F(2,3)*F(3,1)*F(3,2)+F(2,2)^2*F(3,2)*F(3,3)-F(2,2)*F(2,3)*F(3,2)^2);
            f_ratio(2)=-F(3,3)*(F(2,1)*F(3,1)*F(3,3)+F(2,2)*F(3,2)*F(3,3)-F(2,3)*F(3,1)^2-F(2,3)*F(3,2)^2)/ ...
                (F(1,1)*F(1,3)*F(2,1)*F(3,3)-F(1,1)*F(1,3)*F(2,3)*F(3,1)+F(1,2)*F(1,3)*F(2,2)*F(3,3)-F(1,2)*F(1,3)*F(2,3)*F(3,2)+F(2,1)^2*F(2,3)*F(3,3)-F(2,1)*F(2,3)^2*F(3,1)+F(2,2)^2*F(2,3)*F(3,3)-F(2,2)*F(2,3)^2*F(3,2));
            f_ratio=sqrt(f_ratio);
            f_ratio=f_ratio(2)/f_ratio(1);
        case 'Polynomial2'
            f_ratio(1)=-F(3,3)*(F(1,1)*F(1,3)*F(3,3)-F(1,3)^2*F(3,1)+F(2,1)*F(2,3)*F(3,3)-F(2,3)^2*F(3,1))/ ...
                (F(1,1)^2*F(3,1)*F(3,3)+F(1,1)*F(1,2)*F(3,2)*F(3,3)-F(1,1)*F(1,3)*F(3,1)^2-F(1,2)*F(1,3)*F(3,1)*F(3,2)+F(2,1)^2*F(3,1)*F(3,3)+F(2,1)*F(2,2)*F(3,2)*F(3,3)-F(2,1)*F(2,3)*F(3,1)^2-F(2,2)*F(2,3)*F(3,1)*F(3,2));
            f_ratio(2)=-F(3,3)*(F(1,1)*F(3,1)*F(3,3)+F(1,2)*F(3,2)*F(3,3)-F(1,3)*F(3,1)^2-F(1,3)*F(3,2)^2)/ ...
                (F(1,1)^2*F(1,3)*F(3,3)-F(1,1)*F(1,3)^2*F(3,1)+F(1,1)*F(2,1)*F(2,3)*F(3,3)+F(1,2)^2*F(1,3)*F(3,3)-F(1,2)*F(1,3)^2*F(3,2)+F(1,2)*F(2,2)*F(2,3)*F(3,3)-F(1,3)*F(2,1)*F(2,3)*F(3,1)-F(1,3)*F(2,2)*F(2,3)*F(3,2));
            f_ratio=sqrt(f_ratio);
            f_ratio=f_ratio(2)/f_ratio(1);
        case 'Polynomial3'
            % the same as Bougnoux
            f_ratio(1)=-F(3,3)*(F(1,1)*F(1,3)*F(3,2) - F(1,2)*F(1,3)*F(3,1) + F(2,1)*F(2,3)*F(3,2) - F(2,2)*F(2,3)*F(3,1))/(F(1,1)^2*F(3,1)*F(3,2) - F(1,1)*F(1,2)*F(3,1)^2 + F(1,1)*F(1,2)*F(3,2)^2 - F(1,2)^2*F(3,1)*F(3,2) + F(2,1)^2*F(3,1)*F(3,2) - F(2,1)*F(2,2)*F(3,1)^2 + F(2,1)*F(2,2)*F(3,2)^2 - F(2,2)^2*F(3,1)*F(3,2));
            f_ratio(2)=-F(3,3)*(F(1,1)*F(2,3)*F(3,1) + F(1,2)*F(2,3)*F(3,2) - F(1,3)*F(2,1)*F(3,1) - F(1,3)*F(2,2)*F(3,2))/(F(1,1)^2*F(1,3)*F(2,3) - F(1,1)*F(1,3)^2*F(2,1) + F(1,1)*F(2,1)*F(2,3)^2 + F(1,2)^2*F(1,3)*F(2,3) - F(1,2)*F(1,3)^2*F(2,2) + F(1,2)*F(2,2)*F(2,3)^2 - F(1,3)*F(2,1)^2*F(2,3) - F(1,3)*F(2,2)^2*F(2,3));
            f_ratio=sqrt(f_ratio);
            f_ratio=f_ratio(2)/f_ratio(1);
        case 'Polynomial'
            f_ratio = (F(1,1)^2*F(3,1)^2*F(3,3)+2*F(1,1)*F(1,2)*F(3,1)*F(3,2)*F(3,3)-F(1,1)*F(1,3)*F(3,1)^3-F(1,1)*F(1,3)*F(3,1)*F(3,2)^2+F(1,2)^2*F(3,2)^2*F(3,3)-F(1,2)*F(1,3)*F(3,1)^2*F(3,2)-F(1,2)*F(1,3)*F(3,2)^3+F(2,1)^2*F(3,1)^2*F(3,3)+2*F(2,1)*F(2,2)*F(3,1)*F(3,2)*F(3,3)-F(2,1)*F(2,3)*F(3,1)^3-F(2,1)*F(2,3)*F(3,1)*F(3,2)^2+F(2,2)^2*F(3,2)^2*F(3,3)-F(2,2)*F(2,3)*F(3,1)^2*F(3,2)-F(2,2)*F(2,3)*F(3,2)^3)/ ...
                (F(1,1)^2*F(1,3)^2*F(3,3)-F(1,1)*F(1,3)^3*F(3,1)+2*F(1,1)*F(1,3)*F(2,1)*F(2,3)*F(3,3)-F(1,1)*F(1,3)*F(2,3)^2*F(3,1)+F(1,2)^2*F(1,3)^2*F(3,3)-F(1,2)*F(1,3)^3*F(3,2)+2*F(1,2)*F(1,3)*F(2,2)*F(2,3)*F(3,3)-F(1,2)*F(1,3)*F(2,3)^2*F(3,2)-F(1,3)^2*F(2,1)*F(2,3)*F(3,1)-F(1,3)^2*F(2,2)*F(2,3)*F(3,2)+F(2,1)^2*F(2,3)^2*F(3,3)-F(2,1)*F(2,3)^3*F(3,1)+F(2,2)^2*F(2,3)^2*F(3,3)-F(2,2)*F(2,3)^3*F(3,2));
            f_ratio=sqrt(f_ratio);
        case 'Prior'
            return
            % using Prior formula https://cseweb.ucsd.edu/classes/fa04/cse252c/projects/manmohan.pdf
            % TODO - rescale pixel:=0.001
            [f_ratio(1), u1o, v1o, f_ratio(2), u2o, v2o, l1, l2, err, iter] = f1_f2_from_F(F, 2300, 0, 0, 2300, 0, 0, .1, 1, .1, 1);
            f_ratio=f_ratio(2)/f_ratio(1);
        otherwise
            'URGENT: Look at the F2f1f2 function'
    end
else % unit tests
    %p = {[500;600] [500;400]};
    p = {[0;0] [0;0]};
    fg = [1500 900];
    K1 = [fg(1)  0     p{1}(1)
        0      fg(1) p{1}(2)
        0      0     1]
    K2 = [fg(2) 0     p{2}(1)
        0     fg(2) p{2}(2)
        0     0     1]
    % Rotations
    R1 = eye(3)
    R2 = a2r([1;0;0],pi/6)
    % Projection matrices
    P1 = K1*R1*[eye(3)  [-1;0;0]]
    P2 = K2*R2*[eye(3)   [1;0;0]]
    F = (PP2F(P1,P2))
    f_ratio=F2ratio(F,p,'Bougnoux')
    f_ratio=F2ratio(F,p,'Polynomial3')
    f_ratio=F2ratio(F,p)
    fg(2)/fg(1)
    f_ratio = all(abs(f_ratio-fg(2)/fg(1))<1e-8);
    f_ratio=F
end
