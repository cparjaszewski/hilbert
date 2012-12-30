function val = newdocc(fun, c, tol, dck)
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
    % [Krzysztof Parjaszewski, University of Wroclaw]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform used to 
    % better understand and solve the Kramers-Kronig relations in nonlinear optics"
    % krzysztof.parjaszewski@gmail.com

    %% The algorithm:
    
    if nargin < 3, tol = 10^(-3); end;
    if (abs(c) >= 1), error('newdocc:argChk', ['Second argument must be from range (-1,1), but given: ', num2str(c)]); end; 
     
    N = 8;
    lastVal = dosinglecc(fun, c, N/2, dck);
    while (N < 2 ^ 13)
        val = dosinglecc(fun, c, N, dck); 
        N = N .* 2;
        if( abs((val - lastVal) / max(abs([tol, val, lastVal]))) < tol), break; end;
        lastVal = val;
    end
end
    
function val = dosinglecc(fun, c, N, dck)
    l = log2(N);                                   % current iteration count
    d = squeeze(dck(1:N, l))';                     % get d_k coefficients
    
    ks = 0:1:N/2-1; 
    kd = 1 ./ (1 - 4 .* ks .^ 2);                  % denominator in summand value
    ds = d(2 .* ks + 1) .* kd;                     % calculation of d coefficients
    sum1 = sum([0.5, ones(1,N/2-1)] .* ds);        % sumation
    val = 2 .* sum1 + fun(c) .* log((1-c)/(1+c));  % general computation
end