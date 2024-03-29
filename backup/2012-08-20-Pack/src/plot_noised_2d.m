% Plots the CHI3 using extended clenshaw curtis quadrature
function plot_noised_2d(fun,args)
Spectrum = [5.3e-07,5.50e-07,5.70e-07,6.e-07,6.30e-07,6.60e-07, 6.90e-07,7.20e-07,7.50e-07,7.80e-07,8.05e-07,8.50e-07, 8.80e-07,9.20e-07,9.64e-07,1.0150e-06];
Real =     [-3.843580e-36,1.905940e-35,7.855620e-36,3.033850e-37,-1.180530e-36,-3.831240e-36, -6.295560e-36,-5.239590e-36, -9.278530e-36,-1.880560e-35,-5.281310e-36,1.248310e-35, 1.4155e-35,1.0122e-35,5.152260e-36,6.691250e-36];
Imag =     [5.23773E-36, 3.50639E-36, 2.4859E-36, 3.23581E-36, 5.217E-36,8.77457E-36, 1.31922E-35,2.30605E-35,3.47489E-35, 3.12974E-35,1.7376E-35,6.43164E-36,1.13379E-36,9.25422E-38,6.93364E-38,7.19202E-38];

% PEAK-FITTING:
% Gaussian:
% Area = [2.6793503021E-42, 8.997650381479E-43, -1.1774859902702E-41];
% Centers = [7.6375183964162E-7, 6.960970727738E-7, 1.6948356190483E-5];
% Width = [7.0324548928663E-8, 8.7698601291869E-8, 5.2423314696712E-7];
% Height = 2.1783539505151E-36;

%Lorentzian:
Area    = 5.55826E-42;
Centers = 7.55586E-7;
Width   = 9.72123E-8;

% Getting the wavelengths:
minwl = Spectrum(1)/1.8; 
sectrumsize = size(Spectrum);
maxwl = Spectrum(sectrumsize(2))*1.2;

% Prepare arguments for integration
%args = minwl:0.17E-9:maxwl;

% Ranges
gk_a=@(y)(0.99999*y);
gk_b=@(y)(1.00001*y);
a=@(y)(0.9*y);
b=@(y)(1.1*y);

% Lorentzian:
% @(x)(2.*Area/pi).*(Width./(4.*((x-Centers).^2)+Width.^2));  noise
lim=@(x) imag(fun(x)); 
lre=@(x) real(fun(x)); 

% noised: 
xmax = max(args); xmin = min(args); xrange = xmax-xmin;
f=@(x)(1+0.02.*sin(120.*(xmax - x)./xrange)).*fun(x);
g=@(x, y)get_inner_f(f, x, y, a, b);

% Gaussian Sum:
% sum=@(t)0;
% gCount = size(Centers);
% for i=1:gCount(1)
     % sum = @(t)(sum(t) + (Area(i)/(Width(i).*sqrt(pi/2))).*exp(-2.*((t-(Centers(i)))/Width(i)).^2));
% end;

getGK = false;
getTR = true;
getCC = false;
plotRelError = false;

% quadgk
quadgkRefrLower  = @(y)(2/pi) * quadgk(@(x)(x.*f(x)./(x.^2-y.^2)), 0        ,(gk_a(y)));
quadgkRefrUpper  = @(y)(2/pi) * quadgk(@(x)(x.*f(x)./(x.^2-y.^2)),(gk_b(y)) ,      Inf);
quadgkRefrMiddle = @(y)(0);

% quadcc
quadccRefrLower       = @(y) (2/pi) * quadgk(@(x)(x.*f(x)./(x.^2-y.^2)), 0     , (a(y)));
quadccRefrUpper       = @(y) (2/pi) * quadgk(@(x)(x.*f(x)./(x.^2-y.^2)),(b(y)) , Inf);
quadccRefrUpperUpper  = @(y) 0;% (2/pi) * quadgk(@(x)(x.*f(x)./(x.^2-y.^2)),10e3*y ,    Inf);
quadccRefrMiddle      = @(y)(cc_sum(g, y, a(y), b(y)));

tr = @(y)(2/pi)*trapz(args,(innerfun_nlo_g(args,y,f(args))));

% prealocation of values
argsize  = size(args);
valuesCC = 1:argsize(2);
valuesGK = 1:argsize(2);
valuesTR = 1:argsize(2);

% final & main loops

if (getCC==true) 
  for i=1:argsize(2)
    valuesCC(i) = quadccRefrLower(args(i)) + quadccRefrUpper(args(i)) + quadccRefrMiddle(args(i)) + quadccRefrUpperUpper(args(i));
  end;
end;
if (getGK==true) 
  for i=1:argsize(2) 
    valuesGK(i) = quadgkRefrLower(args(i)) + quadgkRefrUpper(args(i)) + quadgkRefrMiddle(args(i));  
  end;
end;
if (getTR==true) 
  for i=1:argsize(2)
    valuesTR(i) = tr(args(i));
  end; 
end;

hold off;
plot(0);
hold on;
if (getGK==true)
  if(plotRelError==true)
    plot(args, (real(fun(args))-valuesGK)./real(fun(args)), 'Color','green');
  else
    plot(args, valuesGK, 'Color','green');
  end;
end;
if (getCC==true)
  if(plotRelError==true)
    plot(args, valuesCC, 'Color','cyan');
  else
    plot(args, valuesCC, 'Color','cyan');
  end;
end;
if (getTR==true)
  if(plotRelError==true)
    plot(args, valuesTR, 'Color','magenta');
  else
    plot(args, valuesTR, 'Color','magenta');
  end;
end;
plot(args, imag(f(args)), 'Color',   'red');
plot(args, real(f(args)), 'Color','yellow');
%plot(Spectrum, Imag, 'Color','blue');

function sum = cc_sum(g, y, a, b)
  % set the n count
  n = 24*2^10;

  % get the chebyshev roots
  cz_sqrs = cos(pi*(0:n)'/n);                                 % Czebyshev roots
  gx = feval(g,cz_sqrs,y)/(2*n);                              % Values of g for roots
  fft_g = real(fft(gx([1:n+1 n:-1:2])));                      % FFT
  c = [fft_g(1); fft_g(2:n)+fft_g(2*n:-1:n+2); fft_g(n+1)];   % Chebyshev indicates
  
  % get the m indicators
  m = 1:1:(n+1);
  d = (2.*y - b - a)/(b-a);
  m(1) = (2/(b-a)) * log(abs((1-d)./(1+d)));
  m(2) = (4/(b-a)) + d * m(1);
  for i=3:1:(n+1)
      if(mod(i,2)==0)
          m(i) = 2*d*m(i-1) - m(i-2);
      else
          m(i) = 2*d*m(i-1) - m(i-2) +  8/((b-a)*(1-(i-1)^2));
      end;
  end;
  
  sum = 1/2 * (c(1) * m(1) + c(2) * m(2));
  for i=3:1:(n+1)
      sum = sum + c(i) * m(i);
  end;
  sum = (b-a)/2 * sum;

function value = get_inner_f(fun,x,y,a,b)
new_x = 0.5*(b(y)-a(y))*x + 0.5*(b(y)+a(y));
value = (new_x) .* fun(new_x)  ./ (new_x + y);