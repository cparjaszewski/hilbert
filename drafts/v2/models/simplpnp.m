function res = simplpnp(delta, Delta)
    % Simple pump-and-probe model
    if nargin < 1, tStart=tic;testhncX2;
        tElapsed=toc(tStart); display(['elapsed time:' num2str(tElapsed)]); 
        finish; 
    end;
    if nargin < 2, Delta = 1.3; end;
    
    G=1;
    gammaba=-0.1; 
    n0=1; 
    eta=1; 
    theta=1.4; 
    Omega1=4.3;
    
    res=G.*n0.*gammaba./(Delta+delta+1i.*eta) .* ...
          (1-Omega1.^2.*(Delta-delta+1i.*eta).*(delta+2.*1i.*eta) ./ ...
          (Delta-1i.*eta)./((delta+1i.*theta).*(Delta+delta+1i.*eta) .* ...
          (delta-Delta+1i.*eta)-Omega1.^2.*(delta+1i.*eta))./2);
end

%% %% TESTS:
%% % A) htran:
%
%% function testhtran1
%     % One-dimensional tests for pnp model (fixed Delta)
%     X = linspace(-30,30,1000);
%     G=1;
%     gammaba=-0.1; 
%     Delta = 1.3;
%     n0=1; 
%     eta=1; 
%     theta=1.4; 
%     Omega1=4.3;
%     
%     Z = G.*n0.*gammaba./(Delta+X+1i.*eta).*(1-Omega1.^2.*(Delta-X+1i*eta).*(X+2.*1i.*eta)./(Delta-1i.*eta)/((X+1i.*theta).*(Delta+X+1i.*eta).*(X-Delta+1i.*eta)-Omega1.^2*(X+1i.*eta))/2);
%     
%     plot(X,real(Z),X,imag(X));
%     plot(X,real(Z),X,imag(Z));
%     plot(X,real(Z),X,imag(Z));
%     G=1;gammaba=-0.1;n0=1;Delta=1.3;eta=1;theta=1.4;Omega1=3;
%     Z = G.*n0.*gammaba./(Delta+X+1i.*eta).*(1-Omega1.^2.*(Delta-X+1i*eta).*(X+2.*1i.*eta)./(Delta-1i.*eta)/((X+1i.*theta).*(Delta+X+1i.*eta).*(X-Delta+1i.*eta)-Omega1.^2*(X+1i.*eta))/2);
%     plot(X,real(Z),X,imag(Z));
%     
%     Z = (G.*n0.*gammaba./(Delta+X+1i.*eta) .* ( 1-Omega1.^2.*(Delta-X+1i.*eta).*(X+2.*1i.*eta)./(Delta-1i.*eta)./((X+1i.*theta).*(Delta+X+1i.*eta).*(X-Delta+1i.*eta)-Omega1.^2.*(X+1i.*eta))./2 ));
%     plot(X,real(Z),X,imag(Z));
%     X = linspace(-10,10,1000);
%     Z = (G.*n0.*gammaba./(Delta+X+1i.*eta) .* ( 1-Omega1.^2.*(Delta-X+1i.*eta).*(X+2.*1i.*eta)./(Delta-1i.*eta)./((X+1i.*theta).*(Delta+X+1i.*eta).*(X-Delta+1i.*eta)-Omega1.^2.*(X+1i.*eta))./2 ));
%     plot(X,real(Z),X,imag(Z));
%     X = linspace(-18,18,1000);
%     Z = (G.*n0.*gammaba./(Delta+X+1i.*eta) .* ( 1-Omega1.^2.*(Delta-X+1i.*eta).*(X+2.*1i.*eta)./(Delta-1i.*eta)./((X+1i.*theta).*(Delta+X+1i.*eta).*(X-Delta+1i.*eta)-Omega1.^2.*(X+1i.*eta))./2 ));
%     [F, HT] = htran(X,Z);
%     plot(X,real(Z),F,HT);
%     % Warning: Imaginary parts of complex X and/or Y arguments ignored 
%     [F, HT] = htran(X,real(Z));
%     plot(X,real(Z),F,HT);
%     plot(X,real(Z),F,HT,X,imag(Z));
%     plot(X,real(Z),F,HT,X,-imag(Z));
%     plot(F,HT,X,-imag(Z));
%     plot(F,HT+imag(Z));
%     plot(F,HT+imag(Z),F,HT,F,-imag(Z));
%     plot(X,real(Z),F,HT+imag(Z),F,HT,F,-imag(Z));
%     clf;
%     plot(X,real(Z),F,HT,F,-imag(Z));
%     hold on, plot(F,HT+imag(Z),'Color','Yellow');
%     clf
%     hold on, plot(F,HT+imag(Z),'Color','Yellow');
%     hold on, plot(F,HT+imag(Z),'Color','Yellow','LineWidth',4);
%     plot(X,real(Z),F,HT,F,-imag(Z));
%     [F1, HT1] = htran(X,imag(Z));
%     hold off;
%     figure;
%     plot(F1,HT1,X,real(Z));
%     plot(X,real(Z),F,HT,F,-imag(Z));
%     plot(X,real(Z),F,HT,F,-imag(Z));
%     hold on, plot(F,HT+imag(Z),'Color','Yellow','LineWidth',4);
%% end
% 
%% function testhtran2
%     % 2-Dimensional tests on the pnp model for htran method
%      X = linspace(-30,30,1000);
%     G=1;
%     gammaba=-0.1; 
%     n0=1; 
%     eta=1; 
%     theta=1.4; 
%     Omega1=4.3;
%     
%     [XX,YY] = meshgrid(X);
%     ZZ = (G.*n0.*gammaba./(YY+XX+1i.*eta) .* ( 1-Omega1.^2.*(YY-XX+1i.*eta).*(XX+2.*1i.*eta)./(YY-1i.*eta)./((XX+1i.*theta).*(YY+XX+1i.*eta).*(XX-YY+1i.*eta)-Omega1.^2.*(XX+1i.*eta))./2 ));
%     
%     mesh(XX,YY,real(ZZ));
%     clf
%     hold off;
%     X = linspace(-9,9,300);
%     [XX,YY] = meshgrid(X);
%     ZZ = (G.*n0.*gammaba./(YY+XX+1i.*eta) .* ( 1-Omega1.^2.*(YY-XX+1i.*eta).*(XX+2.*1i.*eta)./(YY-1i.*eta)./((XX+1i.*theta).*(YY+XX+1i.*eta).*(XX-YY+1i.*eta)-Omega1.^2.*(XX+1i.*eta))./2 ));
%     mesh(XX,YY,real(ZZ));figure,mesh(XX,YY,-imag(ZZ));
%     HTN = zeros(300);
%     for (j=1:100), [~, HTN(j,:)] = htran(XX(j,:),real(ZZ(j,:))); end;
%     mesh(XX,YY,HTN);
%     for (j=1:300), [~, HTN(j,:)] = htran(XX(j,:),real(ZZ(j,:))); end;
%     mesh(XX,YY,HTN);
%     mesh(XX,YY,HTN-imag(ZZ));
%     mesh(XX,YY,HTN+imag(ZZ));
%     for (j=1:300), [~, HTNR(j,:)] = htran(XX(j,:),imag(ZZ(j,:))); end;
%     figure,mesh(XX,YY,HTNR-real(ZZ));
%     mesh(XX,YY,(HTNR-real(ZZ))./real(ZZ));
%     mesh(XX,YY,(HTNR-real(ZZ))./(real(ZZ)+10));
%     mesh(XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+10));
%     mesh(XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+1));
%     mesh(XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+0.1));
%     mesh(XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+0.001));
%     mesh(XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+0.000));
%     mesh(XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+0.01));
%     mesh(XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+2));
%     mesh(XX,YY,(HTNR-real(ZZ)));
%     mesh(XX,YY,(HTNR));
%     mesh(XX,YY,(real(ZZ)));
%     mesh(XX,YY,(imag(ZZ)));
%     mesh(XX,YY,(real(ZZ)));
%     figure,mesh(XX,YY,(HTNR));
%     figure,mesh(XX,YY,(imag(ZZ)));
%     mesh(XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+mean(mean(abs(real(ZZ))))));
%% end
%
%% % B) hncX:
function testhncX1
    % One-dimensional tests for pnp model (fixed Delta)
    noPoints=200;
    ra = -18;
    rb = 18;
    n = 8;
    tol = 10^(-2);
    cs = 0.003;
    
    X = linspace(ra,rb,noPoints);
    Z = simplpnp(X);
    [F,H]=hncX(@(x) real(simplpnp(x)),ra,rb,tol, n, cs,noPoints);   
    figure,plot(X,real(Z),F,-H,F,imag(Z)),title('X,real(Z),F,-H,F,imag(Z),F,-H+imag(Z)');
    hold on, plot(F,-H-imag(Z),'Color','Yellow','LineWidth',4);
    [F1, HT1] = hncX(@(x) imag(simplpnp(x)),ra,rb,tol, n, cs,noPoints); 
    hold off;
    figure,plot(X,imag(Z),F1,HT1,F1,real(Z)),title('X,imag(Z),F1,HT1,F1,real(Z),F1,HT1-real(Z)');
    hold on, plot(F1,HT1-real(Z),'Color','Yellow','LineWidth',4);    
end

function testhncX2
    % 2-Dimensional tests on the pnp model for hncX method
    noPoints=300;
    ra = -30;
    rb = 30;
    n = 8;
    tol = 10^(-2);
    cs = 0.003;

    X = linspace(ra,rb,noPoints);
    [XX,YY] = meshgrid(X);
    ZZ = simplpnp(XX,YY);
    figure,mesh(XX,YY,real(ZZ)),title('XX,YY,real(ZZ)');
    figure,mesh(XX,YY,-imag(ZZ)),title('XX,YY,-imag(ZZ)');
    HTN = zeros(noPoints);
    for j=1:noPoints, [~, HTN(j,:)] = hncX(@(x) real(simplpnp(x,YY(:,j))),ra,rb,tol,n,cs,noPoints);  end;
    figure,mesh(XX,YY,HTN),title('XX,YY,HTN');
    for j=1:noPoints, [~, HTN(j,:)] = hncX(@(x) real(simplpnp(x,YY(:,j))),ra,rb,tol,n,cs,noPoints);  end;
    figure,mesh(XX,YY,HTN),title('XX,YY,HTN');
    figure,mesh(XX,YY,HTN-imag(ZZ)),title('XX,YY,HTN-imag(ZZ)');
    figure,mesh(XX,YY,HTN+imag(ZZ)),title('XX,YY,HTN+imag(ZZ)');
    HTNR=zeros(noPoints);for j=1:noPoints, [~, HTNR(j,:)] =  hncX(@(x) imag(simplpnp(x,YY(:,j))),ra,rb,tol,n,cs,noPoints);  end;
    figure,mesh(XX,YY,HTNR-real(ZZ)),title('XX,YY,HTNR-real(ZZ)');
    figure,mesh(XX,YY,(HTNR-real(ZZ))./real(ZZ)),title('XX,YY,(HTNR-real(ZZ))./real(ZZ)');
    figure,mesh(XX,YY,(HTNR-real(ZZ))./(real(ZZ)+10)),title('XX,YY,(HTNR-real(ZZ))./(real(ZZ)+10)');
    figure,mesh(XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+10)),title('XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+10)');
    figure,mesh(XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+1)),title('XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+1)');
    figure,mesh(XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+0.1)),title('XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+0.1)');
    figure,mesh(XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+0.001)),title('XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+0.001)');
    figure,mesh(XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+0.000)),title('XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+0.000)');
    figure,mesh(XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+0.01)),title('XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+0.01)');
    figure,mesh(XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+2)),title('XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+2)');
    figure,mesh(XX,YY,(HTNR-real(ZZ))),title('XX,YY,(HTNR-real(ZZ))');
    figure,mesh(XX,YY,(HTNR)),title('XX,YY,(HTNR)');
    figure,mesh(XX,YY,(real(ZZ))),title('XX,YY,(real(ZZ))');
    figure,mesh(XX,YY,(imag(ZZ))),title('XX,YY,(imag(ZZ))');
    figure,mesh(XX,YY,(real(ZZ))),title('XX,YY,(real(ZZ))');
    figure,mesh(XX,YY,(HTNR)),title('XX,YY,(HTNR)');
    figure,mesh(XX,YY,(imag(ZZ))),title('XX,YY,(imag(ZZ))');
    figure,mesh(XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+mean(mean(abs(real(ZZ)))))),title('XX,YY,(HTNR-real(ZZ))./(abs(real(ZZ))+mean(mean(abs(real(ZZ))))');
end