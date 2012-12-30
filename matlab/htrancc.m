function H = htrancc(fun, X, tol, inh)
    %% Description:
    % Hilbert transform based on the Clenshaw-Curtis quadrature modified
    % to work on infinite integral which is used in the Hilbert transform formula
        
    %% INPUT:
    % fun - an function_handle                          (function-handle) 
    % X   - interval of abscissas to operate with       (abscissas)
    % tol - deserved tolerance                          (scalar)
    % inh - inner waitbar                               (waitbar-handle)
     
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
    
    a = min(X); b = max(X); NX = length(X); T = abs(b - a);
    A = a - 2 * T; B = b + 2 * T;
    Hh = zeros(1, NX); Y = linspace(a, b, NX);
    
    % Cubic interpolation rutine (N->N-1 points):
    Yp(1)         =   0.375*Y(1)           +   0.75*Y(2)           -  0.125*Y(3);
    Yp(1, NX-1)   =  -0.125*Y(1, NX-2)     +   0.75*Y(1, NX-1)     +  0.375*Y(1, NX);
    Yp(1, 2:NX-2) = -0.0625*Y(1, 1:(NX-3)) + 0.5625*Y(1, 2:(NX-2)) + 0.5625*Y(1, 3:(NX-1)) - 0.0625*Y(1, 4:NX);
    
    % Precalculation of ccFun
    Period = B - A;
    Center = (B+A)/2;
    ccFun = @(t) (fun(t .* Period./2 + Center));
    
    % Preparation of constant values
    AKN = precalculateAKN(ccFun);
    DCK = precalculateDCK(AKN, Yp, A, B);
    
    % Main application loop with the waitbar:
    for n = 1:1:NX-1
        dck = squeeze(DCK(n, :, :));
        Hh(n) = hcc(ccFun, Yp(n), tol, A, B, dck); 
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

function h = hcc(ccFun, y, tol, a, b, dck)
    %% Hilbert transform using the Clenshaw-Curtis quadrature
    
    %% INPUT:
    % ccFun - an function_handle                       (function-handle) 
    % y   - abscissa to operate with                   (scalar)
    % tol - deserved tolerance                         (scalar)
    % a - investigated range minimum                   (scalar)
    % b - investigated range minimum                   (scalar)
    
    %% OUTPUT:
    % h  - the Hilbert transform in the desired point  (scalar)
    
    % We will split the indefinite integral into three main parts:
    % hcc = int(fun(x)/(x-y), x = -inf..a) ...
    %     + int(fun(x)/(x-y), x = a..b) ...
    %     + int(fun(x)/(x-y), x = b..inf)
    % Due to small values ouf of the [a,b] range, we will skip left and
    % right summand.
    
    % Central integral
    h = hccside2side(ccFun, y, tol, a, b, dck);
    
    % Left integral - for now we assume 0
    leftIntVal  = 0;
    
    % Right integral - for now we assume 0
    rightIntVal = 0;
    
    % Final summation
    h = leftIntVal + h + rightIntVal;
end

function h = hccside2side(ccFun, y, tol, a, b, dck)
    %% Hilbert transform using the Clenshaw-Curtis quadrature for a
    %% definite interval
    Period = b - a;
    Center = (b + a) / 2;
    c = (y - Center ) .* 2 ./ Period;
    
    h = 1 / pi * newdocc(@(t) ccFun(t), c, tol, dck);
end

function AKN = precalculateAKN(ccFun)
    % prealocation of memory
    AKN = zeros(2 ^ 13 + 1, 13);
    
    % main loop
    for n=1:13
      N = 2 ^ n;
      akn = precalculateAN(ccFun, N);
      AKN(1:N+1, n) = akn';
    end;
end

function DCK = precalculateDCK(AKN, Yp, a, b)   
    % prealocation of memory
    DCK = zeros(length(Yp), 2 ^ 13 + 2, 13);
    
    Period = b - a;
    Center = (b + a) / 2;
    C = (Yp - Center ) .* 2 ./ Period;
    
    % main loop
    for n=1:13
      N = 2 ^ n;
      a = AKN(1:N+1, n)';
      dck = precalculateDCKN(a, C, N);
      DCK(1:length(C), 1:N+2, n) = dck;
    end;
end

% precalculateAN(fun, N);
% a = get_a(fun, N)
function a = precalculateAN(fun, N) % implementation of (2.1) equation from [1]
  jN  = 0:1:N;                                                        % N+1 points for the last point
  jN1 = 0:1:N-1;                                                      % N points (the last point will be calculated separately)
  args = cos(pi .* jN1 ./ N);                                         % precalculation of arguments
  vals = [0.5, ones(1, N-1), zeros(1,N)] .* [fun(args), zeros(1,N)];  % double prime - first&last element is halved and N additional zeros are added
  ap = 2 / N .* real(fft(vals));                                      % DCT-I done with fft
  al = 1 / N .* fun(-1) .* cos(pi .* jN);                             % FFT has not calculated the last summand
  a = ap(1:N) + al(1:N);                                              % first N values are calculated
  % Last a value:
  argsN = cos(pi .* jN ./ N);                                         % N + 1 range of arguments
  valsN = [0.5, ones(1, N-1), 0.5] .* fun(argsN) .* cos(pi .* jN);    % N + 1 summands
  a(N+1) = 2 / N .*  sum(valsN);                                      % N + 1 th value
end

% precalculateDCKN(A, Yp, N)
% get_d(a, c, N)
function d = precalculateDCKN(a, Yp, N)                             % implementation of (1.11) recursive equation from [1], but vector-style
  l = length(Yp);                                                     % vector length
  d = zeros(l, N+2);                                                  % matrix memory prelocation
  d(:, N+2) = 0; d(:, N+1) = 0;                                       % last and prelast value zeroes
  d(:, N) = a(N+1);                                                   % first one is simple
  for n = N-1:-1:1 
    d(:, n) = 2 * a(n+1) + 2 * Yp' .* d(:, n+1) - d(:, n+2);         % recursive equation solving 
  end; 
end