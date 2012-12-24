function val = newdocc(fun, c, tol)
    %% Description: 
    % Clenshaw-Curtis quadrature to calculate cauchy-principal value
    % integral of the form int(fun(t)/(t-c),-1,1);
    
    %% Based on:
    % [1] T. Hasegawa and T. Torii. An automatic quadrature for Cauchy principal value integrals. Mathematics of Computation, 56:741–754, April 1991.
    
    %% INPUT:
    % fun  -  input function                 (function-handle)
    % c    -  point within the denominator   (scalar)
    % tol  -  required tolerance             (scalar) [default: 10^(-3)]
    
    %% OUTPUT:
    % val  -  integral value (scalar)
    
    %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw, 2011-2012]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform in~nonlinear optics"
    % krzysztof.parjaszewski@gmail.com

    %% The algorithm:
    
    if nargin < 3, tol = 10^(-3); end;
    if (abs(c) >= 1), error('newdocc:argChk', ['Second argument must be from range (-1,1), but given: ', num2str(c)]); end; 
     
    N = 8;
    lastVal = dosinglecc(fun,c,N/2);
    while (N < 1024 * 64)
        val = dosinglecc(fun,c,N); 
        N = N .* 2;
        if( abs((val-lastVal)/max(abs([tol,val,lastVal]))) < tol), break; end;
        lastVal = val;
    end
end
    
function val = dosinglecc(fun, c, N)
    a = get_a(fun, N);                             % get a_k^N coefficients
    d = get_d(a, c, N);                            % get d_k coefficients
     
    ks = 0:1:N/2-1; 
    kd = 1 ./ (1 - 4 .* ks .^ 2);                  % denominator in summand value
    ds = d(2 .* ks + 1) .* kd;                     % calculation of d coefficients
    sum1 = sum([0.5, ones(1,N/2-1)] .* ds);        % sumation
    val = 2 .* sum1 + fun(c) .* log((1-c)/(1+c));  % general computation
end

function a = get_a(fun, N) % implementation of (2.1) equation from [1]
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

function d = get_d(a, c, N)                                               % implementation of (1.11) recursive equation from [1]
  d = zeros(1, N+2);                                                      % memory prelocation
  d(N+2) = 0; d(N+1) = 0;                                                 % last and prelast value zeroes
  d(N) = a(N+1);                                                          % first one is simple
  for n = N-1:-1:1, d(n) = 2 .* a(n+1) + 2 .* c .* d(n+1) - d(n+2);  end; % recursive equation solving
end