% Prepare the 3d plot data for given gaussian approx
function plot_nlo_g_empty()
%% Data for coumarine 307
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
args = minwl:0.17E-9:maxwl;

% Gaussian:
% sum=@(t)0;
% gCount = size(Centers);
% for i=1:gCount(1)
     % sum = @(t)(sum(t) + (Area(i)/(Width(i).*sqrt(pi/2))).*exp(-2.*((t-(Centers(i)))/Width(i)).^2));
% end;

% Lorentzian:
sum=@(t)(2.*Area/pi).*(Width./(4.*((t-Centers).^2)+Width.^2));

refrLower  = @(w)(2/pi) * quadgk(@(t)(t.*sum(t)./(t.^2-w.^2)), 0 ,(w*0.999));
refrUpper  = @(w)(2/pi) * quadgk(@(t)(t.*sum(t)./(t.^2-w.^2)),(w*1.0001) ,Inf);
refrMiddle = @(w)(0); % TODO

tr = @(w)(2/pi)*trapz(args,(innerfun_nlo_g(args,w,sum(args))));

% prealocation of values
argsize  = size(args);
values   = 1:argsize(2);
valuesTR = 1:argsize(2);

% final & main loop
for i=1:argsize(2)
    values(i) = refrLower(args(i)) + refrUpper(args(i)) + refrMiddle(args(i));
    valuesTR(i) = tr(args(i));
end;

hold off;
plot(0,0);
hold on;
plot(args, values,   'Color','green');
plot(args, valuesTR, 'Color','magenta');
plot(args, sum(args),'Color','red');
plot(Spectrum, Real, 'Color','yellow');
plot(Spectrum, Imag, 'Color','blue');