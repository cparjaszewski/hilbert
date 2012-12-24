function [y]=gauss(v,x)
y0=v(1);
A =v(2);
w =v(3);
x0=v(4);
y =y0+(A/w*sqrt(pi/2)).*exp(-2.*(x-x0).^2./w.^2);


