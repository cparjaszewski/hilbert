function H = htrancc(fun, Xinterval, tol, inh)
    %% Description:
    % Hilbert transform based on the Clenshaw-Curtis quadrature modified
    % to work on infinite integral which is used in the Hilbert transform formula
    
    %% Based on:
    % [1] T. Hasegawa and T. Torii. An automatic quadrature for Cauchy principal value integrals. Mathematics of Computation, 56:741–754, April 1991.
      
    %% INPUT:
    % fun       - an function_handle                    (function-handle) 
    % Xinterval - interval of abscissas to operate with (abscissas)
    % tol       - deserved tolerance                    (scalar)
    % inh       - inner waitbar                         (waitbar-handle)
     
    %% OUTPUT:
    % H  - the Hilbert transform on the desired interval (ordinates)
    
    %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform used to 
    % better understand and solve the Kramers-Kronig relations in nonlinear optics"
    % krzysztof.parjaszewski@gmail.com

    %% The algorithm:
    
    % We set the default tolerance if not given:
    if nargin<3, tol=10^(-3); end;
    if nargin<4, h = waitbar(0, 'Please wait', 'Name', 'Hilbert Transform with Clenshaw-Curtis'); else h = inh; end;
    
    a = min(Xinterval); b = max(Xinterval); 
    NX = length(Xinterval); T = abs(b - a);
    A = a - 2 * T; B = b + 2 * T;
    Hh = zeros(1, NX); X = linspace(a, b, NX);
    
    % Cubic interpolation rutine (N->N-1 points):
    Xp(1)         =   0.375*X(1)           +   0.75*X(2)           -  0.125*X(3);
    Xp(1, NX-1)   =  -0.125*X(1, NX-2)     +   0.75*X(1, NX-1)     +  0.375*X(1, NX);
    Xp(1, 2:NX-2) = -0.0625*X(1, 1:(NX-3)) + 0.5625*X(1, 2:(NX-2)) + 0.5625*X(1, 3:(NX-1)) - 0.0625*X(1, 4:NX);
    
    % Precalculation of ccFun
    Period = B - A;
    Center = (B + A) / 2;
    ccFun = @(t) (fun(t .* Period ./2 + Center));
    
    % Preparation of constant values
    AKN = precalculateAKN(ccFun);         %see eq. (1.7) from [1]
    DCK = precalculateDCK(AKN, Xp, A, B); %see eq. (1.11) from [1]
    
    % Main application loop with the waitbar:
    for n = 1:1:NX-1
        dck = squeeze(DCK(n, :, :));
        Hh(n) = hcc(ccFun, Xp(n), tol, A, B, dck); 
        waitbar(n/(NX-1), h, 'Please wait');
    end;
    if nargin<4, close(h); end;
    
    % Reverse cubic interpolation rutine (N-1->N points):
    H = zeros(1,NX);
    H(1)        =  1.875 *Hh(1)        - 1.25  *Hh(2)        + 0.375 *Hh(3);
    H(1,2)      =  0.375 *Hh(1)        + 0.75  *Hh(2)        - 0.125 *Hh(3);
    H(1,3:NX-2) = -0.0625*Hh(1,1:NX-4) + 0.5625*Hh(1,2:NX-3) + 0.5625*Hh(1,3:NX-2) - 0.0625*Hh(1,4:NX-1);
    H(1,NX-1)   = -0.125 *Hh(1,NX-3)   + 0.75  *Hh(1,NX-2)   + 0.375 *Hh(1,NX-1);
    H(1,NX)     =  0.375 *Hh(1,NX-3)   - 1.25  *Hh(1,NX-2)   + 1.875 *Hh(1,NX-1);

end

function h = hcc(ccFun, c, tol, a, b, dck)
    %% Hilbert transform using the Clenshaw-Curtis quadrature
    
    %% INPUT:
    % ccFun - an function_handle                       (function-handle) 
    % c     - abscissa to operate with                 (scalar)
    % tol   - deserved tolerance                       (scalar)
    % a     - investigated range minimum               (scalar)
    % b     - investigated range minimum               (scalar)
    % dck   - precalculated 2d matrix of DCK(n) params (2d matrix)    
    
    %% OUTPUT:
    % h  - the Hilbert transform in the desired point  (scalar)
    
    % We will split the indefinite integral into three main parts:
    % hcc = int(ccFun(x)/(x-c), x = -inf..a) ...
    %     + int(ccFun(x)/(x-c), x = a..b) ...
    %     + int(ccFun(x)/(x-c), x = b..inf)
    % Due to small values ouf of the [a,b] range, we will skip left and
    % right summand.
    
    % Central integral
    h = hccside2side(ccFun, c, tol, a, b, dck);
    
    % Left integral - for now we assume 0
    leftIntVal  = 0;
    
    % Right integral - for now we assume 0
    rightIntVal = 0;
    
    % Final summation
    h = leftIntVal + h + rightIntVal;
end

function h = hccside2side(ccFun, c, tol, a, b, dck)
    %% Description:
    % Hilbert transform using the Clenshaw-Curtis quadrature for a
    % definite interval
    
    %% INPUT:
    % ccFun - an function_handle                       (function-handle) 
    % c     - abscissa to operate with                 (scalar)
    % tol   - deserved tolerance                       (scalar)
    % a     - investigated range minimum               (scalar)
    % b     - investigated range minimum               (scalar)
    % dck   - precalculated 2d matrix of DCK(n) params (2d matrix)    
    
    %% OUTPUT:
    % h  - the Hilbert transform in the desired point  (scalar)
    
    % calculation of cx (based on a, b and c)
    Period = b - a;
    Center = (b + a) / 2;
    cx = (c - Center ) .* 2 ./ Period;
    
    % main calculation
    h = 1 / pi * newdocc(@(t) ccFun(t), cx, tol, dck);
end

function AKN = precalculateAKN(ccFun) % modification of (2.1) equation from [1] to 2d matrix style
    %% Description:
    % Precalculation of AKN params as described in eq. (2.1) from [1]

    %% INPUT:
    % ccFun - an function_handle   (function-handle) 
    
    %% OUTPUT:
    % AKN  - desired matrix        (2d matrix)
    
    % prealocation of memory
    maxIter = getHtranccMaxIter;
    AKN = zeros(2 ^ maxIter + 1, maxIter);
    
    % main loop
    for n=1:maxIter
      N = 2 ^ n;
      akn = precalculateAN(ccFun, N);
      AKN(1:N+1, n) = akn';
    end;
end

function DCK = precalculateDCK(AKN, Xp, a, b) % modification of (1.11) recursive equation from [1] to 3d matrix style
    %% Description:
    % Precalculation of DCK params as described in eq. (1.11) from [1] 

    %% INPUT:
    % AKN   - 2d matrix of AKN coefficients            (2d matrix) [see 1.7 equation from [1]]
    % Xp    - investigated vector of abcissas          (1d vector)
    % a     - investigated range minimum               (scalar)
    % b     - investigated range minimum               (scalar)
    
    %% OUTPUT:
    % DCK  - desired matrix                            (3d matrix)

    % prealocation of memory
    maxIter = getHtranccMaxIter;
    DCK = zeros(length(Xp), 2 ^ maxIter + 2, maxIter);
    
    Period = b - a;
    Center = (b + a) / 2;
    Cp = ((Xp - Center ) .* 2 ./ Period)';
    
    % main loop
    for n=1:maxIter
      N = 2 ^ n;
      a = AKN(1:N+1, n)';
      dck = precalculateDCKN(a, Cp, N);
      DCK(1:length(Cp), 1:N+2, n) = dck;
    end;
end

function a = precalculateAN(ccFun, N) % modification of (2.1) equation from [1] to vector-style
    %% Description:
    % Precalculation of AKN params as described in eq. (2.1) from [1]


    %% INPUT:
    % ccFun - an function_handle        (function-handle) 
    % N     - abscissa to operate with  (scalar)
    
    %% OUTPUT:
    % a  - desired vector               (1d vector)

    jN  = 0:1:N;                                                         % N+1 points for the last point
    jN1 = 0:1:N-1;                                                       % N points (the last point will be calculated separately)
    args = cos(pi .* jN1 ./ N);                                          % precalculation of arguments
    vals = [0.5, ones(1, N-1), zeros(1,N)] .* [ccFun(args), zeros(1,N)]; % double prime - first&last element is halved and N additional zeros are added
    ap = 2 / N .* real(fft(vals));                                       % DCT-I done with fft
    al = 1 / N .* ccFun(-1) .* cos(pi .* jN);                            % FFT has not calculated the last summand
    a = ap(1:N) + al(1:N);                                               % first N values are calculated
    % Last a value:
    argsN = cos(pi .* jN ./ N);                                          % N + 1 range of arguments
    valsN = [0.5, ones(1, N-1), 0.5] .* ccFun(argsN) .* cos(pi .* jN);   % N + 1 summands
    a(N+1) = 2 / N .*  sum(valsN);                                       % N + 1 th value
end

function d = precalculateDCKN(a, Cp, N) % modification of (1.11) recursive equation from [1] to vector-style
    %% Description:
    % Precalculation of DCK params as described in eq. (1.11) from [1] 

    %% INPUT:
    % a  - vector taken from AKN matrix             (1d vector)
    % Cp - abscissas vector scaled to [-1, 1]       (1d vector)
    % N  - number of approximation points           (scalar)
    
    %% OUTPUT:
    % d  - desired matrix                           (2d matrix)

    l = length(Cp);                                                   % vector length
    d = zeros(l, N+2);                                                % matrix memory prelocation
    d(:, N+2) = 0; d(:, N+1) = 0;                                     % last and prelast value zeroes
    d(:, N) = a(N+1);                                                 % first one is simple
    for n = N-1:-1:1 
      d(:, n) = 2 * a(n+1) + 2 * Cp .* d(:, n+1) - d(:, n+2);         % recursive equation solving 
    end; 
end