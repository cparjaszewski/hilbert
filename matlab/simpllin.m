function res = simpllin(omega) 
  % A simple linear function as the Fourier transform of the response
  % signal 1-exp(-t)*sin(20*t)*Theta(t):
  if nargin<1,testhncX;end;
  res= (-20)./(omega.*1i+1-20i)./(omega.*1i+1+20i);
  
end



%% %% TESTS:
%% % A) htran:
%
%% function testhtran
%     X = linspace(-30,30,90);Z = real((-20)./(X.*1i+1-20.*1i)./(X.*1i+1+20.*1i));[~,HT] = htran(X,Z);
%     IZ = -imag(hilbert(Z));PR = -2.236067977.*X.*(-4+5.*X.^2)./(5.*X.^4-9.*X.^2+5);
%     plot(X,HT-PR);
%     figure,plot(X,HT,X,IZ,X,PR);
%     PR = imag((-20)./(X.*1i+1-20.*1i)./(X.*1i+1+20.*1i));
%     plot(X,HT-PR);
%     PR = imag((-20)./(X.*1i+1-20.*1i)./(X.*1i+1+20.*1i));
%     plot(X,HT-PR);
%     plot(X,HT,X,PR);
%     plot(X,HT,X,PR,X,HT-PR);
%     X = linspace(-30,30,1000);Z = real((-20)./(X.*1i+1-20.*1i)./(X.*1i+1+20.*1i));[~,HT] = htran(X,Z);
%     PR = imag((-20)./(X.*1i+1-20.*1i)./(X.*1i+1+20.*1i));plot(X,HT,X,PR,X,HT-PR);
%     X = linspace(-30,30,1000);Z = real((-20)./(X.*1i+1-20.*1i)./(X.*1i+1+20.*1i));[~,HT] = htran(X,Z);
%     PR = imag((-20)./(X.*1i+1-20.*1i)./(X.*1i+1+20.*1i));plot(X,HT,X,PR,X,HT-PR);
%     X = linspace(-30,30,1000);Z = real((-20)./(X.*1i+1-20.*1i)./(X.*1i+1+20.*1i));[~,HT] = htran(X,Z);
%     PR = imag((-20)./(X.*1i+1-20.*1i)./(X.*1i+1+20.*1i));plot(X,HT,X,PR,X,HT-PR);
%     PR = imag((-20)./(X.*1i+1-20.*1i)./(X.*1i+1+20.*1i));plot(X,HT,X,PR,X,HT-PR,'LineWidth',4);
%     hold on,plot(X,Z);
%     clf;
%     hold on;plot(X,HT,X,PR,X,Z);plot(X,HT-PR,'LineWidth',4,'LineColor','Yellow');
%% end
%
%% % B) hncX:
%
function testhncX
    H=hncX(@(x) real(simpllin(x)), -30,30,10^(-2),8,0.003,100);   
    X=linspace(-30,30,100);plot(X,imag(simpllin(X)),X,H);
    H=hncX(@(x) real(simpllin(x)), -30,30,0.2*10^(-1),8,0.03,100);
    plot(X,imag(simpllin(X)),X,H/pi);
    plot(X,imag(simpllin(X)),X,H);
    H=hncX(@(x) real(simpllin(x)), -30,30,0.2*10^(-1),8,0.01,150);
    plot(X,imag(simpllin(X)),X,H);    
    X = linspace(-30,30,150);plot(X,imag(simpllin(X)),X,H);
    X = linspace(-30,30,150);plot(X,real(simpllin(X)),X,H,X,imag(simpllin(X)));hold on; plot(X,H-imag(simpllin(X)),'LineWidth',4,'Color','Yellow');
    figure;
    H=hncX(@(x) -imag(simpllin(x)), -30,30,0.2*10^(-1),8,0.01,150);
    X = linspace(-30,30,150);plot(X,imag(simpllin(X)),X,H,X,real(simpllin(X)));hold on; plot(X,H-real(simpllin(X)),'LineWidth',4,'Color','Yellow');
    figure,H=hncX(@(x) real(simpllin(x)), -30,30,0.2*10^(-1),3,0.01,150);plot(X,real(simpllin(X)),X,H,X,imag(simpllin(X)));hold on; plot(X,H-imag(simpllin(X)),'LineWidth',4,'Color','Yellow');
    figure,H=hncX(@(x) -imag(simpllin(x)), -30,30,0.2*10^(-1),3,0.01,150);plot(X,imag(simpllin(X)),X,H,X,real(simpllin(X)));hold on; plot(X,H-real(simpllin(X)),'LineWidth',4,'Color','Yellow');
    figure,H=hncX(@(x) real(simpllin(x)), -30,30,0.2*10^(-1),13,0.5,150);plot(X,real(simpllin(X)),X,H,X,imag(simpllin(X)));hold on; plot(X,H-imag(simpllin(X)),'LineWidth',4,'Color','Yellow');
    H=hncX(@(x) real(simpllin(x)), -30,30,0.2*10^(-1),13,0.05,150);plot(X,real(simpllin(X)),X,H,X,imag(simpllin(X)));hold on; plot(X,H-imag(simpllin(X)),'LineWidth',4,'Color','Yellow');
    figure,H=hncX(@(x) -imag(simpllin(x)),-30,30,0.2*10^(-1),13,0.05,150);plot(X,imag(simpllin(X)),X,H,X,real(simpllin(X)));hold on; plot(X,H-real(simpllin(X)),'LineWidth',4,'Color','Yellow');
    figure,H=hncX(@(x) real(simpllin(x)), -30,30,0.2*10^(-2),13,0.05,150);plot(X,real(simpllin(X)),X,H,X,imag(simpllin(X)));hold on; plot(X,H-imag(simpllin(X)),'LineWidth',4,'Color','Yellow');
    figure,H=hncX(@(x) -imag(simpllin(x)),-30,30,0.2*10^(-1),23,0.05,150);plot(X,imag(simpllin(X)),X,H,X,real(simpllin(X)));hold on; plot(X,H-real(simpllin(X)),'LineWidth',4,'Color','Yellow');
    figure,H=hncX(@(x) real(simpllin(x)), -30,30,0.2*10^(-2),23,0.05,150);plot(X,real(simpllin(X)),X,H,X,imag(simpllin(X)));hold on; plot(X,H-imag(simpllin(X)),'LineWidth',4,'Color','Yellow');
    figure,H=hncX(@(x) -imag(simpllin(x)),-30,30,0.2*10^(-1),33,0.05,150);plot(X,imag(simpllin(X)),X,H,X,real(simpllin(X)));hold on; plot(X,H-real(simpllin(X)),'LineWidth',4,'Color','Yellow');
    figure,H=hncX(@(x) real(simpllin(x)), -30,30,0.2*10^(-2),33,0.05,150);plot(X,real(simpllin(X)),X,H,X,imag(simpllin(X)));hold on; plot(X,H-imag(simpllin(X)),'LineWidth',4,'Color','Yellow');

end