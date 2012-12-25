% Gauss
figure;
t = (-15:1/5:30);
x1 = gaussmf(t,[2 5])+0.04*(randn(size(t))-0.5);
x2 = gaussmf(t,[2 5]);
y1 = hilbert(x1);
y2 = hilbert(x2);
plot(t,real(y1),'.','Color','red'), hold on;
plot(t,real(y2),'Color','red');
plot(t,imag(y2));
plot(t,imag(y1),':'), title('Gauss'), hold off;

% Lorentzian
figure;
x1 = 1./(1+(t-3).^2)+0.04*(randn(size(t))-0.5);
x2 = 1./(1+(t-3).^2);
y1 = hilbert(x1);
y2 = hilbert(x2);
plot(t,real(y1),'.','Color','red'), hold on;
plot(t,real(y2),'Color','red');
plot(t,imag(y2));
plot(t,imag(y1),':'), title('Lorentz'), hold off;


% nonlinear K-K
figure;
omega_O = 1; T_2 = 1; s = 1; A = 1;
x1 = A .* ((t - omega_O) + 1i./T_2)./(t.^2 - 2i.*t./T_2.*(1+s.^2).^0.5 + (omega_O.^2+1./T_2.^2.*(1+s.^2)));
% x2 = A * ((t - omega_O) + 1i./T_2)/(t.^2 - 2i.*t./T_2.*(1+s^2).^0.5 + (omega_O^2+1/T_2.^2.*(1+s.^2)))+0.004*(randn(size(t))-0.5);
y1 = hilbert(x1);
% y2 = hilbert(x2);
plot(t,real(y1),'.','Color','red'), hold on;
%plot(t,real(y2),'Color','red');
%plot(t,imag(y2));
plot(t,imag(y1),':'), title('nonlinear K-K'), hold off;
