
   
    [F,A,support(i)]=F_features(u1,u2,'|F|=0',testset,3,false);
    F=reshape(F,3,3);
    f(i,:)=F2f1f2(F);
    f(i,:)=f(i,:)*diag([1/A{1}(1) 1/A{2}(1)]).*norm_(i,:);
    
    %vanilla
    % pixel in this case = 1000, therefore i scale down the error so that
    % it in the same scale as that of nister-solved one.
    original(i)=sampson_error(F,testset{1},testset{2})/1000;
    %original(i)=sampson_error(F,u1,u2);
    
    % corrected
    [corrF,K1,K2]=correctF(F,abs(f(i,:)));
    corrected(i)=sampson_error(corrF,testset{1},testset{2})/1000;
    if norm(F/F(1)-corrF/corrF(1),'fro')>1e-6 & (abs(imag(f(i,1)))<eps) & (abs(imag(f(i,2)))<eps)
        norm(F/F(1)-corrF/corrF(1),'fro')
        f(i,:)
        norm_(i,:)
        A{1}
        A{2}
        warning('calcFocals');
    end
    % recomputed
    u1t=inv(K1)*a2h(u1);
    u2t=inv(K2)*a2h(u2);
    testt={inv(K1)*a2h(testset{1}) inv(K2)*a2h(testset{2})};
    Es=E5ptNister([(u1t(:,1:5)); (u2t(:,1:5))]);
    Es=cell2mat(reshape(Es,1,1,[]));
    recomputed(i)=sampson_error(voteF(Es,nan,testt,0.001),testt{1},testt{2});
    
    % method used by Tomas in frecexample
    PP = E2PP(K2'*F*K1);
    E  = PP2F(PP{1}{1},PP{2}{1});
    weird(i)=sampson_error(E,testt{1},testt{2});
  