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
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform in~nonlinear optics"
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
    
    % Main application loop with the waitbar:
    for n = 1:1:NX-1
        Hh(n) = hcc(fun, Yp(n), tol, A, B); 
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

function h = hcc(fun, y, tol, a, b)
    %% Hilbert transform using the Clenshaw-Curtis quadrature
    
    %% INPUT:
    % fun - an function_handle                         (function-handle) 
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
    h = hccside2side(fun, y, tol, a, b);
    
    % Left integral
    leftIntVal  = 0;
    
    % Right integral
    rightIntVal = 0;
    
    % Final summation
    h = leftIntVal + h + rightIntVal;
end

function h = hccside2side(fun, y, tol, a, b)
    %% Hilbert transform using the Clenshaw-Curtis quadrature for a
    %% definite interval
    Period = b - a;
    Center = (b+a)/2;
    d = (y - Center ) .* 2 ./ Period;
   
    ccFun = @(t) (fun(t .* Period./2 + Center));
    
    h = 1 / pi * newdocc(@(t) ccFun(t), d, tol);
end