clear;
x  = linspace(-2,2,400);
y = x;
often = 60;
intensity = 0.4;
[XX,YY] = meshgrid(x,y);
fun = @(r,j) (exp(-1i.*(r+1i.*j)).*(1-intensity.*sin(often.*r)).*(1+intensity.*cos(often.*j))); Z = fun(XX,YY);
figure;
subplot(1,2,1);mesh(XX,YY,imag(Z),real(Z));subplot(1,2,2);mesh(XX,YY,real(Z),imag(Z));

figure;