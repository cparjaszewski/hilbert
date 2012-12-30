function res = quadncX(fun,a,b,tol,n)
    %% Description:
    % A simple algorithm of Newton-Cotes is presented.
    % We do not like to use Newton-Cotes of much higher degree, instead we use the
    % complex quadrature of lower degree for many small steps
    
    %% INPUT:
    % fun - the function we will like to itegrate (function-handle)
    % a - a integration starting point
    % b - integration ending point
    % tol - required tollerance
    % n - degree of the Newton-Cotes quadrature
    
    %% OUTPUT:
    % res  - calculated quadrature value
    
    %% Author info:
		% [Krzysztof Parjaszewski, University of Wroclaw]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform in~nonlinear optics"
    % krzysztof.parjaszewski@gmail.com
    
    %% The algorithm:
    
    if nargin < 3, error('QUADNCX:ParamErr', 'Wrong number of parameter given: fun, a ,b are mandatory'); end;
    if nargin < 4, tol = 10^(-5); end;
    if nargin < 5, n = 8; end;
    
    % While this method does not shows very good accuraccy, a set of
    % arbitrary set parameters has been prepared, but they may be changed
    % in case of different quadratures.
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
    
%     if (iterates >= maxIterate), warning('QUADNCX:MaxIterReached', 'maximum number of iterations reached - the integral seems to be singular');
%     elseif (nosteps >=maxNoSteps), warning('QUADNCX:MaNoSteps', 'maximum number of steps reached - the integral seems to be singular'); 
%     end;
end

function tab = prevaluateNCTab(n,a,b)
    % The table of NC coefficient is calculated using Maple:
    % Digits:=100;n:=(n);evalf(seq(1/n*(-1)^(n-k)*1/(factorial(k)*factorial(n-k))*int(product(t-j,j=0..n)/(t-k),t=0..n),k=0..n));
    maplestr = ['Digits:=40:n:=' num2str(n) ':tab:=evalf(seq(((' num2str(b) ')-(' num2str(a) '))/n*(-1)^(n-k)*1/(factorial(k)*factorial(n-k))*int(product(t-jj,jj=0..n)/(t-k),t=0..n),k=0..n));'];
    maple string;
    tab = str2num(maple(maplestr))';
end

function res = doStep(fun,a,b,nosteps,tab)
    % A simple numerical quadrature for already calculated coefficients if
    % performed
    nc = size(tab,1);
    intervals = linspace(a,b,nosteps);
    
    % Quadrature is divided into small steps
    M = zeros(nc,nosteps-1);
    M(1,:) = intervals(1:nosteps-1);
    M(nc,:) = intervals(2:nosteps);
    h = (M(nc,:) - M(1,:))./nc;
    
    % Final calculation in two lines
    for k=2:(nc-1),  M(k,:) = M(1,:) + (k-1).*h; end;
    res = sum(fun(M') * tab)/nosteps;
end