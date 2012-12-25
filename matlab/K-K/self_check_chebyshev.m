function self_check_chebyshev()
%% TODO!

A=3.1;
T1=0.4;
fun = @(x) (1 ./ (A - x - 1i.*T1));
%args=0.01:0.003:8;
args=0.02:0.02:4;
[X,Y] = meshgrid(args,args);

a=@(y)(0.9.*y);
b=@(y)(1.3.*y);
f=@(x)imag(fun(x)); 

new_x = @(x,y) (0.5.* x .* (b(y) - a(y)) + 0.5.*(b(y) + a(y)));
jedna = @(x,y) ((new_x(x, y) .* f(new_x(x, y))) ./ (new_x(x, y) + y));

Zjedna = jedna(X,Y);
druga = @(x,y) (drugaFun(x,y,jedna));
argsize = size(args);
Zdruga = X;
for k=1:1:argsize(2)
    for j=1:1:argsize(2)
%       Zdruga(k,j) = druga(X(j),Y(k));
    end;
end;
%Z = abs(Zjedna - Zdruga)./Zjedna;
hold off;
plot(0);
%hold on;
mesh(X, Y, Zjedna);
  
function value = drugaFun(x,y,g)
% set the n count
n = 48;

% get the chebyshev roots
cz_sqrs = cos(pi*(0:n)'/n);                                 % Czebyshev roots
gx = feval(g,cz_sqrs,y)/(2*n);                              % Values of g for roots
fft_g = real(fft(gx([1:n+1 n:-1:2])));                      % FFT
c = [fft_g(1); fft_g(2:n)+fft_g(2*n:-1:n+2); fft_g(n+1)];   % Chebyshev indicates

T = 1:1:(n+1);
T(1) = 1;
T(2) = x;
for i=3:1:(n+1)
  Tprev = T(i-1);
  Tprevprev = T(i-2);
  T(i) = (2.*x.*Tprev - Tprevprev);
end;

Tone = T(1);
Tlast = T(n+1);
sum = 1/2 * (c(1) * Tone); % first coefficient is divided by two
for i=2:1:n
  sum = sum + c(i) * T(i);
end;
sum = sum + 1/2 * (c(n+1) * Tlast); % last coefficient is divided by two

% return value
value = sum;