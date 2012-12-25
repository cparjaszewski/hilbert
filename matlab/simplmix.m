function res = simplmix(delta,Delta)
    if nargin<2,Delta = 1; end;

    N=1;
    w0=1;
    muba=1;
    T1=1;
    T2=2;
    Omega2=1;
    epsilon0=1;
    h=1;

    res = 2./3.*N.*w0.*muba.^4.*(-delta-Delta-1i./T2).*(delta+2.*1i./T2) ./ ...
          (Delta+1i./T2)./epsilon0./(h.^3)./(Delta+delta+1i./T2) ./ ...
             (delta.^3 - ...
              2.*1i.*delta.^2./T2 - ...
              delta.*Delta.^2 - ...
              2.*1i./T1.*delta./T2 - ...
              delta./(T2.^2) + ...
              delta.^2./T1 - ...
              Delta.^2./T1-1./(T1.*T2.^2) - ...
              Omega2.*delta + ...
              Omega2./T2.*1i);

end

%% %% TESTS
%% % a) htran
%
%% function testhtran1
%     % One-dimensional tests for mix model (fixed Delta)
%     X = linspace(-9,9,300);
%     Delta = 1; 
%     N=1;
%     w0=1;
%     muba=1;
%     T1=1;
%     T2=2;
%     Omega2=1;
%     epsilon0=1;
%     h=1;
%     
%     Z = 2./3.*N.*w0.*muba.^4.*(-X-Delta-1i./T2).*(X+2.*1i./T2)./(Delta+1i./T2)./epsilon0./(h.^3)./(Delta+X+1i./T2)./(X.^3-2.*1i.*X.^2./T2-X.*Delta.^2-2.*1i./T1.*X./T2-X./(T2.^2)+X.^2./T1-Delta.^2./T1-1./(T1.*T2.^2)-Omega2.*X+Omega2./T2.*1i);
%     [F,HT] = htran(X,real(Z));
%     plot(X,real(Z),F,HT);
%     plot(X,real(Z),F,HT,X,imag(Z));
%     hold on;plot(HT-imag(Z),'Color','Yellow','LineWidth',4);
%     hold off;
%     plot(X,real(Z),F,HT,X,imag(Z));
%     hold on;plot(X,HT-imag(Z),'Color','Yellow','LineWidth',4);
%     X = linspace(-9,9,1000);
%     Z = 2./3.*N.*w0.*muba.^4.*(-X-Delta-1i./T2).*(X+2.*1i./T2)./(Delta+1i./T2)./epsilon0./(h.^3)./(Delta+X+1i./T2)./(X.^3-2.*1i.*X.^2./T2-X.*Delta.^2-2.*1i./T1.*X./T2-X./(T2.^2)+X.^2./T1-Delta.^2./T1-1./(T1.*T2.^2)-Omega2.*X+Omega2./T2.*1i);
%      [F,HT] = htran(X,real(Z));
%     hold off;
%     plot(X,real(Z),F,HT,X,imag(Z));
%     hold on;plot(X,HT-imag(Z),'Color','Yellow','LineWidth',4);
%% end

%% function testhtran2
%     % 2-Dimensional tests on the mix model for htran method
%     N=1;
%     w0=1;
%     muba=1;
%     T1=1;
%     T2=2;
%     Omega2=1;
%     epsilon0=1;
%     h=1;
%     
%     X = linspace(-9,9,300);
%     [XX,YY] = meshgrid(X);
%     ZZ = 2./3.*N.*w0.*muba.^4.*(-XX-YY-1i./T2).*(XX+2.*1i./T2)./(YY+1i./T2)./epsilon0./(h.^3)./(YY+XX+1i./T2)./(XX.^3-2.*1i.*XX.^2./T2-XX.*YY.^2-2.*1i./T1.*XX./T2-XX./(T2.^2)+XX.^2./T1-YY.^2./T1-1./(T1.*T2.^2)-Omega2.*XX+Omega2./T2.*1i);
%     mesh(XX,YY,real(ZZ));
%     hold off;
%     mesh(XX,YY,real(ZZ));
%     figure,mesh(XX,YY,imag(ZZ));
%     MTN=zeros(300);for j=1:300, [~, MTN(j,:)] = htran(XX(j,:),real(ZZ(j,:))); end;
%     MTNR=zeros(300);for j=1:300, [~, MTNR(j,:)] = htran(XX(j,:),imag(ZZ(j,:))); end;
%     figure, mesh(XX,YY,MTN);
%     mesh(XX,YY,MTN-real(ZZ));
%     mesh(XX,YY,MTN-imag(ZZ));
%     figure,mesh(XX,YY,MTNR-real(ZZ));
%     figure,mesh(XX,YY,MTNR+real(ZZ));
%     mesh(XX,YY,(MTNR-real(ZZ))/(abs(real(ZZ))+mean(mean(real(ZZ)))));
%     mesh(XX,YY,(MTNR-real(ZZ))./(abs(real(ZZ))+mean(mean(real(ZZ)))));
%     mesh(XX,YY,(MTNR-real(ZZ))./(abs(real(ZZ)+mean(mean(real(ZZ))))));
%     mesh(XX,YY,(MTNR-real(ZZ))./(abs(real(ZZ)+mean(mean(abs(real(ZZ)))))));
%     mesh(XX,YY,MTNR+real(ZZ));
%     mesh(XX,YY,(MTNR+real(ZZ))./((abs(real(ZZ)+mean(mean(abs(real(ZZ))))))));
%     mesh(XX,YY,((abs(real(ZZ)+mean(mean(abs(real(ZZ))))))));
%     mesh(XX,YY,(MTNR+real(ZZ)));
%     figure,mesh(XX,YY,(MTNR));
%     mesh(XX,YY,(MTNR+real(ZZ))./((mean(mean(abs(real(ZZ)))))));
%     mesh(XX,YY,(MTNR+real(ZZ))./(abs(real(ZZ))+(mean(mean(abs(real(ZZ)))))));
%     mesh(XX,YY,(MTNR+real(ZZ))./(abs(real(ZZ)+(mean(mean(abs(real(ZZ))))))));
%     mesh(XX,YY,(MTNR+real(ZZ))./(abs(real(ZZ))+(mean(mean(abs(real(ZZ)))))));
%     MTNR2=zeros(300,300);for j=1:300, [~, MTNR2(j,:)] = htran(XX(j,:),imag(ZZ(j,:))); end;
%     mesh(XX,YY,XX);
%     mesh(XX,YY,YY);
%% end