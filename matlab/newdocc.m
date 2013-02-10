function val = newdocc(ccFun, cx, tol, dck)
    %% Description: 
    % Clenshaw-Curtis quadrature to calculate cauchy-principal value
    % integral of the form int(ccFun(t)/(t-cx),-1,1);
    
    %% Based on:
    % [1] T. Hasegawa and T. Torii. An automatic quadrature for Cauchy principal value integrals. Mathematics of Computation, 56:741–754, April 1991.
    
    %% INPUT:
    % ccFun - input function, scaled to [-1, 1] (function-handle)
    % cx    - point within the denominator      (scalar)
    % dck   - precalculated mx of cc params     (2d matrix)
    % tol   - required tolerance                (scalar) [default: 10^(-3)]
    
    %% OUTPUT:
    % val  -  integral value (scalar)
    
    %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform in~nonlinear optics"
    % krzysztof.parjaszewski@gmail.com

    %% The algorithm:
    
    if nargin < 4, tol = 10^(-3); end;
    if (abs(cx) >= 1), error('newdocc:argChk', ['Second argument must be from range (-1,1), but given: ', num2str(cx)]); end; 
     
    N = 8;
    lastVal = dosinglecc(ccFun, cx, N/2, dck);
    maxIter = getHtranccMaxIter;
    while (N < 2 ^ maxIter)
        val = dosinglecc(ccFun, cx, N, dck); 
        N = N .* 2;
        if( abs((val - lastVal) / max(abs([tol, val, lastVal]))) < tol), break; end;
        lastVal = val;
    end
end
    
function val = dosinglecc(ccFun, cx, N, dck)
    l = log2(N);                                      % current iteration count, see eq. (1.12) from [1]
    d = squeeze(dck(1:N, l))';                        % get d_k coefficients, see eq. (1.11) from [1]
    
    ks = 0:1:N/2-1; 
    kd = 1 ./ (1 - 4 .* ks .^ 2);                     % denominator in summand value
    ds = d(2 .* ks + 1) .* kd;                        % calculation of d coefficients
    sum1 = sum([0.5, ones(1,N/2-1)] .* ds);           % sumation
    val = 2 .* sum1 + ccFun(cx) .* log((1-cx)/(1+cx));  % general computation, see eq. (1.10) from [1]
end