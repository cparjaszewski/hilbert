function res=pertdata2(omega1, omega2, conf)
    %% Description:
    % Source file used for modeling the simple, second-order quantum-perturbative model
    
    %% INPUT:
    % omega1 - discrete array of X values  (vector)
    % omega2 - discrete array of Y values  (vector)
    % conf   - model configuration         (special)
    
    %% OUTPUT:
    % res - comple array of model values   (vector)

    %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform used to 
    % better understand and solve the Kramers-Kronig relations in nonlinear optics"
    % krzysztof.parjaszewski@gmail.com

    %% Test case:
    % if no arguments are given, tests are run
    if (nargin==0)
        tic
        testhncX2; 
        res = toc;
        msg = ['Second-order quantum-perturbative model test took: ' num2str(res) ' seconds'];
        display(msg);
    else
        N     = conf.N;     % for example:  5;
        e0    = conf.e0;    % for example:  1;
        h     = conf.h;     % for example: -1;
        M     = conf.M;     % for example: [[1,2];[-1,-2]];
        Omega = conf.Omega; % for example: [[3, 16];[4, 12]];
        gamma = conf.gamma; % for example: [[1,2];[-1,3]];

        s = 0;
        for l=1:2,
          for m=1:2,
              for n=1:2,
                s = s + M(l,n).*M(n,m).*M(m,l)./(Omega(n,l)-omega1-omega2-1i.*gamma(n,l))./(Omega(m,l)-omega1-1i.*gamma(m,l)) ...
                      + M(l,n).*M(n,m).*M(m,l)./(Omega(n,l)-omega1-omega1-1i.*gamma(n,l))./(Omega(m,l)-omega2-1i.*gamma(m,l)) ...
                      + M(l,n).*M(n,m).*M(m,l)./(Omega(m,n)-omega1-omega2-1i.*gamma(m,n))./(Omega(n,l)+omega2+1i.*gamma(n,l)) ...
                      + M(l,n).*M(n,m).*M(m,l)./(Omega(m,n)-omega1-omega2-1i.*gamma(m,n))./(Omega(n,l)+omega2+1i.*gamma(n,l)) ... 
                      + M(l,n).*M(n,m).*M(m,l)./(Omega(n,m)+omega1+omega2+1i.*gamma(n,m))./(Omega(m,l)-omega1-1i.*gamma(m,l)) ...
                      + M(l,n).*M(n,m).*M(m,l)./(Omega(n,m)+omega1+omega2+1i.*gamma(n,m))./(Omega(m,l)-omega1-1i.*gamma(m,l)) ...
                      + M(l,n).*M(n,m).*M(m,l)./(Omega(m,l)+omega1+omega2+1i.*gamma(m,l))./(Omega(n,l)+omega1+1i.*gamma(n,l)) ...
                      + M(l,n).*M(n,m).*M(m,l)./(Omega(m,l)+omega1+omega2+1i.*gamma(m,l))./(Omega(n,l)+omega2+1i.*gamma(n,l));
             end;
          end;
        end;
        res = s .* N./2./e0./(h.^2);
    end
end
%% %% TESTS:

function testhtran
    % 2-Dimensional test for htran method for the pertdata2 model
    conf.M     = [[1,2];[-1,-2]];
    conf.Omega = [[3, 16];[4, 12]];
    conf.gamma = [[1,2];[-1,3]];
    X = linspace(-30,30,500); [XX,YY] = meshgrid(X); ZZ = pertdata2(XX, YY, conf);
    clf,mesh(XX,YY,real(ZZ)),figure,mesh(XX,YY,imag(ZZ));
    QP2R=zeros(500);for j=1:500 ,QP2R(j,:) = htran(XX(j,:),real(ZZ(j,:))); end;
    QP2I=zeros(500);for j=1:500 ,QP2I(j,:) = htran(XX(j,:),imag(ZZ(j,:))); end;
    figure,mesh(XX,YY,QP2I),figure,mesh(XX,YY,QP2R);
    for j=1:500 ,QP2I(:,j) = htran(XX(:,j),real(ZZ(:,j))); end;
    for j=1:500 ,QP2I(:,j) = htran(XX(:,j),imag(ZZ(:,j))); end;
    for j=1:500 ,QP2I(:,j) = htran(YY(:,j),imag(ZZ(:,j))); end;
    for j=1:500 ,QP2R(:,j) = htran(YY(:,j),imag(ZZ(:,j))); end;
    for j=1:500 ,QP2I(:,j) = htran(YY(:,j),real(ZZ(:,j))); end;
    figure,mesh(XX,YY,QP2I),figure,mesh(XX,YY,QP2R);
    for j=1:500 ,[~,QP2R(:,j)] = htran(YY(:,j),imag(ZZ(:,j))); end;
    for j=1:500 ,[~,QP2I(:,j)] = htran(YY(:,j),real(ZZ(:,j))); end;
    figure,mesh(XX,YY,QP2I),figure,mesh(XX,YY,QP2R);
    for j=1:500 ,[~,QP2R(j,:)] = htran(XX(j,:),imag(ZZ(j,:))); end;
    for j=1:500 ,[~,QP2I(j,:)] = htran(XX(j,:),real(ZZ(j,:))); end;
    figure,mesh(XX,YY,QP2I),figure,mesh(XX,YY,QP2R);
    figure,mesh(XX,YY,QP2I),title('IKK'),figure,mesh(XX,YY,QP2R),title('RKK'),figure,mesh(XX,YY,real(ZZ)),title('real'),figure,mesh(XX,YY,imag(ZZ)),title('imag');
    figure,mesh(XX,YY,-QP2I),title('IKK'),figure,mesh(XX,YY,QP2R),title('RKK'),figure,mesh(XX,YY,real(ZZ)),title('real'),figure,mesh(XX,YY,imag(ZZ)),title('imag');
    figure,mesh(XX,YY,QP2I),title('IKK'),figure,mesh(XX,YY,-QP2R),title('RKK'),figure,mesh(XX,YY,real(ZZ)),title('real'),figure,mesh(XX,YY,imag(ZZ)),title('imag');
end

function testhncX2
    %% Description:
    % 2-Dimensional tests on the perturbative 2-D model for hncX method
    
    %% Quadrature configuration:
    pts = 60;
    a = -30;
    b = 30;
    n = 8;
    tol = 10^(-1);
    cs = 0.3;
    wrn = false; % we don't want the warnings here
    
    %% Model configuration:
    conf.N  = 8;
    conf.e0 = 1.4;
    conf.h  = -2.7;
    conf.M  = [[3, -0.5]; [1.2, 2.4]];
    conf.Omega  = [[-3, 13]; [5, 9.7]];
    conf.gamma  = [[0.7, 2.3]; [-4, -1.3]];

    %% Main calculations:
    X = linspace(a, b, pts);
    [XX,YY] = meshgrid(X); %% input 2-D grid of points
    ZZ = pertdata2(XX, YY, conf); % calculated model 
    figure,mesh(XX,YY,real(ZZ)),title('XX,YY,real(ZZ)');
    figure,mesh(XX,YY,-imag(ZZ)),title('XX,YY,-imag(ZZ)');
    HNCX2 = zeros(pts);
    for j=1:pts, [~, HNCX2(j,:)] = hncX(@(x) real(pertdata2(x, YY(:,j), conf)),a,b,tol,n,cs,pts, wrn);  end;
    figure,mesh(XX,YY,HTN),title('XX,YY,HNCX2');
    for j=1:pts, [~, HNCX2(j,:)] = hncX(@(x) real(pertdata2(x, YY(:,j), conf)),a,b,tol,n,cs,pts, wrn);  end;
    figure,mesh(XX,YY,HNCX2),title('XX,YY,HNCX2');
    figure,mesh(XX,YY,HNCX2-imag(ZZ)),title('XX,YY,HNCX2-imag(ZZ)');
    figure,mesh(XX,YY,HNCX2+imag(ZZ)),title('XX,YY,HNCX2+imag(ZZ)');
    HNCRX2=zeros(pts);
    for j=1:pts, [~, HNCRX2(j,:)] =  hncX(@(x) imag(pertdata2(x, YY(:,j), conf)),a,b,tol,n,cs,pts, wrn);  end;
    figure,mesh(XX,YY,HNCRX2-real(ZZ)),title('XX,YY,HNCRX2-real(ZZ)');
    figure,mesh(XX,YY,(HNCRX2-real(ZZ))./real(ZZ)),title('XX,YY,(HNCRX2-real(ZZ))./real(ZZ)');
    figure,mesh(XX,YY,(HNCRX2-real(ZZ))./(real(ZZ)+10)),title('XX,YY,(HNCRX2-real(ZZ))./(real(ZZ)+10)');
    figure,mesh(XX,YY,(HNCRX2-real(ZZ))./(abs(real(ZZ))+10)),title('XX,YY,(HNCRX2-real(ZZ))./(abs(real(ZZ))+10)');
    figure,mesh(XX,YY,(HNCRX2-real(ZZ))./(abs(real(ZZ))+1)),title('XX,YY,(HNCRX2-real(ZZ))./(abs(real(ZZ))+1)');
    figure,mesh(XX,YY,(HNCRX2-real(ZZ))./(abs(real(ZZ))+0.1)),title('XX,YY,(HNCRX2-real(ZZ))./(abs(real(ZZ))+0.1)');
    figure,mesh(XX,YY,(HNCRX2-real(ZZ))./(abs(real(ZZ))+0.001)),title('XX,YY,(HNCRX2-real(ZZ))./(abs(real(ZZ))+0.001)');
    figure,mesh(XX,YY,(HNCRX2-real(ZZ))./(abs(real(ZZ))+0.000)),title('XX,YY,(HNCRX2-real(ZZ))./(abs(real(ZZ))+0.000)');
    figure,mesh(XX,YY,(HNCRX2-real(ZZ))./(abs(real(ZZ))+0.01)),title('XX,YY,(HNCRX2-real(ZZ))./(abs(real(ZZ))+0.01)');
    figure,mesh(XX,YY,(HNCRX2-real(ZZ))./(abs(real(ZZ))+2)),title('XX,YY,(HNCRX2-real(ZZ))./(abs(real(ZZ))+2)');
    figure,mesh(XX,YY,(HNCRX2-real(ZZ))),title('XX,YY,(HNCRX2-real(ZZ))');
    figure,mesh(XX,YY,(HNCRX2)),title('XX,YY,(HNCRX2)');
    figure,mesh(XX,YY,(real(ZZ))),title('XX,YY,(real(ZZ))');
    figure,mesh(XX,YY,(imag(ZZ))),title('XX,YY,(imag(ZZ))');
    figure,mesh(XX,YY,(real(ZZ))),title('XX,YY,(real(ZZ))');
    figure,mesh(XX,YY,(HNCRX2)),title('XX,YY,(HNCRX2)');
    figure,mesh(XX,YY,(imag(ZZ))),title('XX,YY,(imag(ZZ))');
    figure,mesh(XX,YY,(HNCRX2-real(ZZ))./(abs(real(ZZ))+mean(mean(abs(real(ZZ)))))),title('XX,YY,(HNCRX2-real(ZZ))./(abs(real(ZZ))+mean(mean(abs(real(ZZ))))');
end