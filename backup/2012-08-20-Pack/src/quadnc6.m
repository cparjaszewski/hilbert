function res = quadnc6(fun,a,b,tol)
    % We do not use Newton-Cotes of much higher degree, instead we use the
    % complex quadrature of lower degree for many small steps
    if nargin < 3, error('QUADSN:ParamErr', 'Wrong number of parameter given: fun, a ,b are mandatory'); end;
    if nargin < 4, tol = 100000*eps; end;
    nosteps = 100; maxNoSteps = 1000000;
    iterates = 1; maxIterate = 10000;
    M = 1.1;
    while (iterates < maxIterate);
        nStep = ceil(nosteps * M);
        step1 = doStep(fun,a,b,nosteps);
        step2 = doStep(fun,a,b,nStep);
        res = (4*step2-step1)/3; 
        err = abs(step1-res)/1000;
        if (err<tol), break;
        else nosteps = nStep; M = M*1.3;
        end;
        if (nosteps >=maxNoSteps), break; end;
    end
    if (iterates >= maxIterate), warning('QUADSN:MaxIterReached', 'maximum number of iterations reached - the integral seems to be singular');
    elseif (nosteps >=maxNoSteps), warning('QUADSN:MaNoSteps', 'maximum number of steps reached - the integral seems to be singular'); 
    end;
end



function res = doStep(fun,a,b,nosteps)
    % This table is for NC-6 calculated within Maple:
    % Digits:=100;n:=6;evalf(seq(1/n*(-1)^(n-k)*1/(factorial(k)*factorial(n
    % -k))*int(product(t-j,j=0..n)/(t-k),t=0..n),k=0..n));
    tab = [0.4880952380952380952380952380952380952380952380952380952380952380952380952380952380952380952380952381e-1, ...
           0.2571428571428571428571428571428571428571428571428571428571428571428571428571428571428571428571428571, ...
           0.3214285714285714285714285714285714285714285714285714285714285714285714285714285714285714285714285714e-1, ...
           0.3238095238095238095238095238095238095238095238095238095238095238095238095238095238095238095238095238, ...
           0.3214285714285714285714285714285714285714285714285714285714285714285714285714285714285714285714285714e-1, ...
           0.2571428571428571428571428571428571428571428571428571428571428571428571428571428571428571428571428571, ... 
           0.4880952380952380952380952380952380952380952380952380952380952380952380952380952380952380952380952381e-1]';
       
    nc = size(tab,1);
    h = (b-a)/nosteps;
    intervals = linspace(a,b,nosteps);
    
    % Simpson rule for small steps
    M = zeros(nc,nosteps-1);
    M(1,:) = intervals(1:nosteps-1);
    M(nc,:) = intervals(2:nosteps);
    h = (M(nc,:) - M(1,:))./nc;
    
    for i=2:(nc-1),  M(i,:) = M(1,:) + (i-1).*h;    end;
    res = sum(fun(M') * tab)/nosteps;
end