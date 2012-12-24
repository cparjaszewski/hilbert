% Plots the CHI3
function plot_chi_3()
  clear all;
  mSize   = 100;
  x       = linspace(0.3,5,mSize);
  y       = x;
  [XX,YY] = meshgrid(x,y);
  ZRef3   = XX;
  ZAbs3   = XX;
  ZRef3T  = XX;
  ZAbs3T  = XX;
  A       = 2.1;
  B       = 0.0;
  C       = 2.9;
  T1      = 0.2; 
  T2      = -0.3;
  T3      = 0.1;
  Chi31   = @(w1,w2) (1 ./ ((A - w1 - 1i.*T1).*(B           - 1i.*T2).*(C - w2 - 1i.*T3)));
  Chi32   = @(w1,w2) (1 ./ ((A + w1 + 1i.*T1).*(B           - 1i.*T2).*(C - w2 - 1i.*T3)));
  Chi33   = @(w1,w2) (1 ./ ((A + w1 + 1i.*T1).*(B + w1 - w2 + 1i.*T2).*(C - w2 - 1i.*T3)));
  Chi34   = @(w1,w2) (1 ./ ((A + w1 + 1i.*T1).*(B + w1 - w2 + 1i.*T2).*(C + w1 + 1i.*T3)));
  Chi3    = @(w1,w2) Chi31(w1,w2)+Chi32(w1,w2)+Chi33(w1,w2)+Chi34(w1,w2);
  Abs3    = @(w1,w2) imag(Chi3(w1,w2));
  Ref3    = @(w1,w2) real(Chi3(w1,w2));
  Z       = Chi33(XX, YY);
  near    = 0.001;

  CalcRef3Lower = @(w1,w2) ((2./pi)      .* quadgk(@(t) t .* Abs3(t,w2) ./ (t.^2-w1.^2) , 0, w1*(1-near)));
  CalcRef3Upper = @(w1,w2) ((2./pi)      .* quadgk(@(t) t .* Abs3(t,w2) ./ (t.^2-w1.^2) ,    w1*(1+near), Inf));
  CalcAbs3Lower = @(w1,w2) ((-2.*w1./pi) .* quadgk(@(t)      Ref3(t,w2) ./ (t.^2-w1.^2) , 0, w1*(1-near)));
  CalcAbs3Upper = @(w1,w2) ((-2.*w1./pi) .* quadgk(@(t)      Ref3(t,w2) ./ (t.^2-w1.^2) ,    w1*(1+near), Inf));
  % Trapez method
  CalcRef3Trapz = @(w1,w2) ((2/pi)       .* trapz(x, (x .* innerfun_chi_g(Abs3(x,w2),x,w1))));
  CalcAbs3Trapz = @(w1,w2) ((-2.*w1./pi) .* trapz(x, (     innerfun_chi_g(Ref3(x,w2),x,w1))));


h = waitbar(0,'Please Wait');
for i = 1:mSize
    for j = 1:mSize
        % Numerically evaluate integral, adaptive Gauss-Kronrod quadrature:
         ZRef3(j,i) = CalcRef3Lower(x(i),y(j)) + CalcRef3Upper(x(i),y(j));
         ZAbs3(j,i) = CalcAbs3Lower(x(i),y(j)) + CalcAbs3Upper(x(i),y(j));
        
        % Trapezoidal integration:
         ZAbs3T(j,i) = CalcAbs3Trapz(x(i),y(j));
         ZRef3T(j,i) = CalcRef3Trapz(x(i),y(j));
    end;
    waitbar(i/mSize,h);
end;
close(h);

tr = figure;
title('Trapez');

qu = figure;
title('Quad');

% Values - th. refraction , colors - th. absorption
figure(qu);
subplot(2,3,1);
mesh(XX,YY,real(Z),imag(Z));
title('QU raw values');

figure(tr);
subplot(2,3,1);
mesh(XX,YY,real(Z),imag(Z));
title('TR raw values');

% Relative error for real part:
m = max(max(real(Z)));

figure(qu);
subplot(2,3,2);
mesh(XX,YY,abs(ZRef3-real(Z))./m);
title('QU absolute error - real part');

figure(tr);
subplot(2,3,2);
mesh(XX,YY,abs(ZRef3T-real(Z))./m);
title('TR absoluteerror - real part');

% Relative error for real part:
m = max(max(imag(Z)));
figure(qu);
subplot(2,3,3);
mesh(XX,YY,abs(ZAbs3-imag(Z))./m);
title('QU absolute error - imag part');

figure(tr);
subplot(2,3,3);
mesh(XX,YY,abs(ZAbs3T-imag(Z))./m);
title('TR absolute error - imag part');

% K-K calculated Absorption
figure(qu);
subplot(2,3,4);
mesh(XX,YY,ZAbs3);
title('QU calculated imag part');

figure(tr);
subplot(2,3,4);
mesh(XX,YY,ZAbs3T);
title('TR calculated imag part');

% K-K calculated Refraction
figure(qu);
subplot(2,3,5);
mesh(XX,YY,ZRef3);
title('QU calculated real part');

figure(tr);
subplot(2,3,5);
mesh(XX,YY,ZRef3T);
title('TR calculated real part');

figure(qu);
subplot(2,3,6);
mesh(XX,YY,abs(((ZAbs3-imag(Z))./ZAbs3)),abs(ZAbs3-imag(Z)));
title('QU relative/absolute error(v/c)- imag');

figure(tr);
subplot(2,3,6);
mesh(XX,YY,abs(((ZAbs3T-imag(Z))./ZAbs3T)),abs(ZAbs3T-imag(Z)));
title('TR relative/absolute error(v/c)- imag');


% Theoretical absorption plot

figure;
subplot(1,2,1);
mesh(XX,YY,min(5*ones(size(ZAbs3)),abs((ZAbs3T - ZAbs3)./ZAbs3)));
subplot(1,2,2);
mesh(XX,YY,min(5*ones(size(ZAbs3)),abs((ZRef3T - ZRef3)./ZRef3)));
%mesh(XX,YY,imag(Z));
%mesh(XX,YY,real(Z),imag(Z));