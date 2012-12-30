function [HT,F] = htran(F,R)
    %% Description:
    % Hilbert transform for discrete value vector. This is a fast and quite satisfactory 
    % algorithm based on the important assumption that R disappears out of the F region.

    %% Based on:
    % 1. "Hilbert Transform by Numerical Integration" by I.J. Weinbert, 
    % RADC-TR-79-3, January 1979, Rome Aire Development Center
    % ERRATA:
    % In the original publication there was an typo error - the Hilbert
    % transform is not the 
    %   1/pi*int(u(tau)/(tau-t),tau=-inf..inf) 
    % but should rather be:
    %   1/pi*int(u(tau)/(t-tau),tau=-inf..inf) 
    
       
    %% INPUT:
    % F - an array of points (abscissas)
    % R - an array of values (ordinates)
    
    %% OUTPUT:
    % HT   - an array of Hilbert Transform values * (ordinates)
    % F    - an array of points *                   (abscissas) 
    % *    - corrected, so size(F,1) <= size(F,2)
    
    % We also do not assume that the input signal R is even or odd - because
    % the original Hilbert Transform was prepared for all type of possible
    % functions.

    %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform in~nonlinear optics"
    % krzysztof.parjaszewski@gmail.com

    %% The algorithm:
   
    % Correcing the input size:
    if(size(R,1)>size(R,2)), R=R';end;
    if(size(F,1)>size(F,2)), F=F';end;
    N = size(R,2);
    
    % Calculating the F middlepoint array:
    Fp(1,1:N-1) = 0.5 .* (F(1,1:N-1)+F(1,2:N));
    step = F(2)-F(1);
    
    % Build h array (step array for faster calculation) - Simpson's rule combined with the trapezoidal rule
    h(1) = step/3;
    h(1,2:N-1) = 2.*step./3 + mod(1:N-2,2).*2.*step./3;
    h(N) = step/3;
    if (mod(N,2)==1), h(N-1) = 5*step/6; h(N) = step/2; end;
    
    % Cubic interpolation for the interior values and parabolic interpolation for the end values we obtain
    % new set of ordinates, to omit the singularity:
    Rp = zeros(1,N-1);
    Rp(1)       =   0.375*R(1)      +  0.75*R(2)      - 0.125*R(3);
    Rp(1,N-1)   =  -0.125*R(1,N-2)  +  0.75*R(1,N-1)  + 0.375*R(1,N);
    Rp(1,2:N-2) = -0.0625*R(1,1:N-3)+0.5625*R(1,2:N-2)+0.5625*R(1,3:N-1)-0.0625*R(1,4:N);
    
    % Build Yp array - the heart and first step in integration:
    Yp = zeros(1,N-1);
    for i=1:N-1, Yp(i) = sum(h.*(R-Rp(i))./(F-Fp(i))); end;
    
    % Build Xp array - the second step in integration:
    Xp = -Rp./pi .* log((Fp-Fp(1))./(Fp(N-1)-Fp)) + Yp./pi;
    
    % Build Xres array - translating the Xp array into Xres array (cubic
    % interpolation, as previous) - but this time Xp is smaller than Xres:
    Xres = zeros(1,N);
    Xres(1)       =  1.875 *Xp(1)       - 1.25  *Xp(2)       + 0.375 *Xp(3);
    Xres(1,2)     =  0.375 *Xp(1)       + 0.75  *Xp(2)       - 0.125 *Xp(3);
    Xres(1,3:N-2) = -0.0625*Xp(1,1:N-4) + 0.5625*Xp(1,2:N-3) + 0.5625*Xp(1,3:N-2) - 0.0625*Xp(1,4:N-1);
    Xres(1,N-1)   = -0.125 *Xp(1,N-3)   + 0.75  *Xp(1,N-2)   + 0.375 *Xp(1,N-1);
    Xres(1,N)     =  0.375 *Xp(1,N-3)   - 1.25  *Xp(1,N-2)   + 1.875 *Xp(1,N-1);
    
    % We take the real/imag separately:
    % Orig = imag(Xres);
    HT = Xres;
end

