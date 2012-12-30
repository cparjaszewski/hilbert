function [F, H] = hncX(fun, a, b, tol, n, cs, pts, wrn)
    %% Description:
    % This algorithm is a simple implementation of the Hilbert transform based on the Newton-Cotes quadrature
    % of an arbitrary degree with getting close, but just omitting the singularity
    % We do not like to use Newton-Cotes of much higher degree, instead we use the
    % complex quadrature of lower degree for many small steps  
    
    %% INPUT:
    % fun   - real-valued function to transform    (function-handle)
    % a     - integration starting point           (scalar)
    % b     - integration ending point             (scalar)
    % tol   - required tolerance                   (scalar)  [optional, default: 10^(-5)]
    % n     - degree of Newton-Cotes quadrature    (scalar)  [optional, default: 8]
    % cs    - how close to singularity             (scalar)  [optional, default: 0.01]
    % pts   - #points to perform Hilbert transform (scalar)  [optional, default: 200]
    % wrn   - warning on flag                      (boolean) [optional, defalt: true]
    
    %% OUTPUT:
    % F - an array of points                       (abscissas) 
    % H - calculated Hilbert transform values      (ordinates)
    
    %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform used to 
    % better understand and solve the Kramers-Kronig relations in nonlinear optics"
    % krzysztof.parjaszewski@gmail.com
    
    %% The algorithm:
    
    % Preparation of arguments
    if nargin < 3, error('QUADNCX:ParamErr', 'Wrong number of parameter given: fun, a, b are mandatory'); end;
    if nargin < 4, tol = 10^(-5); end;
    if nargin < 5, n = 8; end;
    if nargin < 6, cs = 0.01; end;
    if nargin < 7, pts = 200; end;
    if nargin < 8, wrn = true; end;
    
    % We would like to integrate within 5 times bigger region
    aMax = a - 4*abs(a); bMax = b+4*abs(b);
    F = linspace(a,b,pts);
    H = zeros(1,pts);
    
    % Not to make user wait too long, we are showing the waitbar
    hBar=waitbar(0,'Please Wait');
    
    % The main hilbert transform loop
    for k=1:pts
        waitbar(k/pts);
        innerfun = @(x) fun(x)./(x-(F(k)));
        H(k) = quadncX(@(x)innerfun(x),aMax,F(k)-cs,tol,n, wrn) + quadncX(@(x)innerfun(x),F(k)+cs,bMax,tol,n, wrn);
    end
    % finally we are closing the waitbar
    close(hBar);
    
    H = H/pi;  
end

function res = quadncX(fun, a, b, tol, n, wrn)
    %% Description:
    % Newton-Cotes quadrature of an arbitrary degree

    %% More info:
    % While this method does not shows very good accuraccy, a set of
    % arbitrary set parameters has been prepared, but they may be changed
    % in case of different quadratures.
    
    %% INPUT:
    % fun   - real-valued function to integrate     (function-handle)
    % a     - integration starting point            (scalar)
    % b     - integration ending point              (scalar)
    % tol   - required tolerance                    (scalar)  [optional, default: 10^(-5)]
    % n     - degree of the Newton-Cotes quadrature (scalar)  [optional, default: 8]
    % wrn   - warning on flag                       (boolean) [optional, defalt: true]
    
    %% OUTPUT:
    % res - an quadrature result point              (scalar)
    
    nosteps = 100; maxNoSteps = 100000;
    iterates = 1; maxIterate = 10000;
    M = 1.1;
    tab = prevaluateNCTab(n,a,b);
    
    while (iterates < maxIterate);
        nStep = ceil(nosteps * M);
        step1 = doStep(fun,a,b,nosteps,tab);
        step2 = doStep(fun,a,b,nStep,tab);
        res = (4*step2-step1)/3; 
        err = abs(step1-res);
        if (err<tol), break;
        else nosteps = nStep; M = M*1.3;
        end;
        if (nosteps >=maxNoSteps), break; end;
    end
    if (wrn == true)
     if (iterates >= maxIterate), warning('QUADNCX:MaxIterReached', 'maximum number of iterations reached - the integral seems to be singular');
     elseif (nosteps >=maxNoSteps), warning('QUADNCX:MaNoSteps', 'maximum number of steps reached - the integral seems to be singular'); 
     end;
    end
end

function tab = prevaluateNCTab(n,a,b)
    %% Description:
    % A routine to pre-evalute the Newton-Cotes coefficients

    %% More info:
    % The table of NC coefficient is calculated using Maple syntax:
    % Digits:=100;n:=(n);evalf(seq(1/n*(-1)^(n-k)*1/(factorial(k)*factorial(n-k))*int(product(t-j,j=0..n)/(t-k),t=0..n),k=0..n));
    
    %% INPUT:
    % n     - degree of the Newton-Cotes quadrature (scalar) 
    % a     - integration starting point            (scalar)
    % b     - integration ending point              (scalar)
    
    %% OUTPUT:
    % tab - Newton-Cotes coefficients               (vector)
    
    %% The algorithm:
    maplestr = ['Digits:=40:n:=' num2str(n) ':tab:=evalf(seq(((' num2str(b) ')-(' num2str(a) ...
        '))/n*(-1)^(n-k)*1/(factorial(k)*factorial(n-k))*int(product(t-jj,jj=0..n)/(t-k),t=0..n),k=0..n));'];
    maple string;
    tab = str2num(maple(maplestr))'; %#ok<ST2NM>
end

function res = doStep(fun, a, b, nosteps, tab)
    %% Description:
    % Simple numerical quadrature for given coefficients is performed
    
    %% INPUT:
    % fun     - the real-valued function we will like to integrate (function-handle)
    % a       - integration starting point                         (scalar)
    % b       - integration ending point                           (scalar)
    % nosteps - number of points to divine the integration range   (scalar)
    % tab     - vector of Newton-Cotes coeffients                  (vector) 
    
    %% OUTPUT:
    % res     - calculated quadrature value                        (scalar)

    %% The algorithm:
    nc = size(tab,1);
    intervals = linspace(a, b, nosteps);
    
    % Quadrature is divided into small steps:
    M = zeros(nc, nosteps-1);
    M(1,:) = intervals(1:nosteps-1);
    M(nc,:) = intervals(2:nosteps);
    h = (M(nc,:) - M(1,:))./nc;
    
    % Final calculation in two lines but using matrix-multiplication:
    for k=2:(nc-1),  M(k,:) = M(1,:) + (k-1).*h; end;
    res = sum(fun(M') * tab)/nosteps;
end