function H = htrancc(fun,X,tol)
    %% Description:
    % Hilbert transform based on the Clenshaw-Curtis quadrature modified
    % to work on infinite integral which is used in the Hilbert transform formula
        
    %% INPUT:
    % fun - an function_handle                          (function-handle) 
    % X   - interval of abscissas to operate with       (abscissas)
    % tol - deserved tolerance                          (scalar)
     
    %% OUTPUT:
    % H  - the Hilbert transform on the desired interval (ordinates)
    
    %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw, Summer 2011]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform used to 
    % better understand and solve the Kramers-Kronig relations in nonlinear optics"
    % krzysztof.parjaszewski@gmail.com

    %% The algorithm:
    
    % We set the default tolerance if not given:
    if nargin<3,tol=10^(-3); end;
    
    a = min(X);b=max(X); NX = length(X); 
    Hh = zeros(1,NX); Y = linspace(a,b,NX);
    
    % Cubic interpolation rutine (N->N-1 points):
    Yp(1)        =   0.375*Y(1)       +  0.75*Y(2)       - 0.125*Y(3);
    Yp(1,NX-1)   =  -0.125*Y(1,NX-2)  +  0.75*Y(1,NX-1)  + 0.375*Y(1,NX);
    Yp(1,2:NX-2) = -0.0625*Y(1,1:(NX-3))+0.5625*Y(1,2:(NX-2))+0.5625*Y(1,3:(NX-1))-0.0625*Y(1,4:NX);
    
    % Main application loop with the waitbar:
    h = waitbar(0,'Please wait','Name','Hilbert Transform with Clenshaw-Curtis');
    for n = 1:1:NX-1
        Hh(n) = hcc(fun,Yp(n), tol); 
        waitbar(n/(NX-1),h,'Please wait');
    end;
    close(h);
    
    % Reverse cubic interpolation rutine (N-1->N points):
    H = zeros(1,NX);
    H(1)        =  1.875 *Hh(1)        - 1.25  *Hh(2)        + 0.375 *Hh(3);
    H(1,2)      =  0.375 *Hh(1)        + 0.75  *Hh(2)        - 0.125 *Hh(3);
    H(1,3:NX-2) = -0.0625*Hh(1,1:NX-4) + 0.5625*Hh(1,2:NX-3) + 0.5625*Hh(1,3:NX-2) - 0.0625*Hh(1,4:NX-1);
    H(1,NX-1)   = -0.125 *Hh(1,NX-3)   + 0.75  *Hh(1,NX-2)   + 0.375 *Hh(1,NX-1);
    H(1,NX)     =  0.375 *Hh(1,NX-3)   - 1.25  *Hh(1,NX-2)   + 1.875 *Hh(1,NX-1);

end

function h = hcc(fun, y, tol)
    %% Hilbert transform using the Clenshaw-Curtis quadrature
    
    %% INPUT:
    % fun - an function_handle                         (function-handle) 
    % y   - abscissa to operate with                   (scalar)
    % tol - deserved tolerance                         (scalar)
    
    %% OUTPUT:
    % h  - the Hilbert transform in the desired point  (scalar)
    
    % We need to change variable in the infinite integration x =
    % t/(1-t.^2) and we assume limit(fun(t)/t)=0 for +-infinity as should
    % be in Hilbert transform
    ccFun = @(t) (abs(t)<1.*fun(t./(1-t.^2)).*(1+t.^2)./((t.^2-1).*(t+y.*t.^2-y)));
    
    % We set N to the initial value:
    N = 8; lastVal = 0;
    while (true)
        h = 1/pi*docc(ccFun,N);
        
        % A simple adaptive method - if the relative error is to large - we
        % double the points count:
        if (abs(lastVal-h)/(abs(h)+tol)<tol),break; end;
        N = N*2;
        lastVal = h;
    end
end
