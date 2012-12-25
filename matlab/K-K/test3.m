function test3
    % Gaussian
    figure;
    t1 = [9.85222E-4 0.00104 0.00109 0.00114 0.00118 0.00124 0.00129 0.00133 0.00139 0.00145 0.00152 0.00159 0.00167 0.00175 0.00182 0.00189];
    r1 = [-3.84358E-36 1.90594E-35 7.85562E-36 3.03385E-37 -1.18053E-36 -3.83124E-36 -6.29556E-36 -5.23959E-36 -9.27853E-36 -1.88056E-35 -5.28131E-36 1.24831E-35 1.4155E-35 1.0122E-35 5.15226E-36 6.69125E-36];
    gaussPlot(t1,r1,t1);

    figure;
    t2 = [9.85222E-4:0.0000104:0.00189];
    gaussPlot(t2,r1,t1);
    
    % Lorentzian
    
end


function gaussPlot(t,r,t1)
    A = [4.46828E-39 2.36367E-39];
    t0 = [0.00131 0.00144];
    w = [1.22998E-4	2.21996E-4];
    xa1 = A(1)./(w(1).*sqrt(pi./2)).*exp(-2.*(t - t0(1)).^2./w(1).^2);
    xa2 = A(2)./(w(2).*sqrt(pi./2)).*exp(-2.*(t - t0(2)).^2./w(2).^2);
    x1 = xa1+xa2;
    x2 = x1 + 1.0E-35.*0.1.*(randn(size(t))-0.5);
    y1 = hilbert(x1);
    y2 = hilbert(x2);
    plot(t,real(y1),'.','Color','red'), hold on;
    plot(t,real(y2),'Color','red');
    plot(t1,r,'Color','green');
    plot(t,imag(y2));
    plot(t,imag(y1),':'), title('Gauss'), hold off;
end



% Lorentzian
% figure;
% x1 = 1./(1+(t-3).^2)+0.04*(randn(size(t))-0.5);
% x2 = 1./(1+(t-3).^2);
% y1 = hilbert(x1);
% y2 = hilbert(x2);
% plot(t,real(y1),'.','Color','red'), hold on;
% plot(t,real(y2),'Color','red');
% plot(t,imag(y2));
% plot(t,imag(y1),':'), title('Lorentz'), hold off;
% 
% 
% % nonlinear K-K
% figure;
% omega_O = 1; T_2 = 1; s = 1; A = 1;
% x1 = A .* ((t - omega_O) + 1i./T_2)./(t.^2 - 2i.*t./T_2.*(1+s.^2).^0.5 + (omega_O.^2+1./T_2.^2.*(1+s.^2)));
% % x2 = A * ((t - omega_O) + 1i./T_2)/(t.^2 - 2i.*t./T_2.*(1+s^2).^0.5 + (omega_O^2+1/T_2.^2.*(1+s.^2)))+0.004*(randn(size(t))-0.5);
% y1 = hilbert(x1);
% % y2 = hilbert(x2);
% plot(t,real(y1),'.','Color','red'), hold on;
% %plot(t,real(y2),'Color','red');
% %plot(t,imag(y2));
% plot(t,imag(y1),':'), title('nonlinear K-K'), hold off;
