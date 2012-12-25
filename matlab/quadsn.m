function res = quadsn(fun, a, b, tol)
    % We do not use Newton-Cotes of much higher degree, instead we use the
    % complex quadrature of lower degree (for smaller steps)
    if nargin < 3, error('QUADSN:ParamErr', 'Wrong number of parameter given: fun, a ,b are mandatory'); end;
    if nargin < 4, tol = 100000*eps; end;
    nosteps = 100; maxNoSteps = 100000000;
    iterates = 1; maxIterate = 10000;
    M = 1.1;
    while (iterates < maxIterate);
        nStep = ceil(nosteps * M);
        step1 = doStep(fun,a,b,nosteps);
        step2 = doStep(fun,a,b,nStep);
        err = abs(step1-step2)/1000;
        if (err<tol)
            res = (4*step2-step1)/3; 
        break;
        else nosteps = nStep;
        end;
        if (nosteps >=maxNoSteps), break; end;
    end
    if (iterates >= maxIterate), warning('QUADSN:MaxIterReached', 'maximum number of iterations reached - the integral seems to be singular');
    elseif (nosteps >=maxNoSteps), warning('QUADSN:MaNoSteps', 'maximum number of steps reached - the integral seems to be singular'); end;
end

function res = doStep(fun,a,b,nosteps)
    h = (b-a)/nosteps;
    intervals = a:h:b;
    
    % Simpson rule for small steps
    first = intervals(1:nosteps-1); val1 = fun(first);
    second = intervals(2:nosteps); val2 = fun(second);
    third = (first+second)./2; val3 = 4.*fun(third);
    
    % sum it all
    res = sum(val1 + val2 + val3)/(6*nosteps); 
end