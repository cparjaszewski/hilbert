function res = hiltrans(fun, s, tol)
    % We would like to make a hilbert transform of a weighted function H(wf;x)
    % = int(f(t)/(t-x)*w(t),t=-infinity..infinity);
    % 
    % 1. Find points x(m,w) with which interpolate function f into matrix
    % F(m,w)
    % 2. Calculate matrix A(m,w) 
    % 3. Calculate A.*F

    % Howto 2. ?
    % 2.a. Calculate matrix lambda(m,w) with Christoffel constants with respect to the weight w,
    % 2.b. Calculate polynomial p_j into matrix P(j,m,w)
    % 2.c. Calculate polynomial q_j into matrix Q(j) a little bit redundantly
    % H(wp_j(w))

    % Howto obtain more convergency and the stability?
    % -.1562500000e-1*(-.2195602390e1507*s*exp(4.*s^2+235.6194490*I*s)+14.17963
    % 

    if(nargin<3), tol = 10^(-6); end;
    newFun = @(t) fun(t)./(t-s);
    left = getLeftHT(newFun, s, tol);
    middle = getMiddleHT(newFun, s, tol);
    right = getRightHT(newFun, s, tol);
    
    res = 1./pi * (left + middle + right);
end

function left = getLeftHT(fun, s, tol)
    alpha = s - 25;
    beta = s-0.1;
    
    iterates = 1; maxIterate = 400;
    sum = 0; prev = quadnc6(fun, alpha, beta, tol);
    while (iterates< maxIterate), iterates = iterates+1;
        sum = sum + prev;
        
        pAlpha = alpha;
        alpha = alpha - 2.*(beta-alpha);
        beta = pAlpha;
        
        next = quadnc6(fun, alpha, beta, tol);
        if ((next/prev < tol) || (abs(prev-next)<tol)), break; end;
        prev = next;
    end;
    
    % Something gone wrong
    if (iterates >= maxIterate), warning('HILTRANSgetLeftHT:MaxIterReached', 'maximum number of iterations reached - the integral seems to be singular'); end;
    
    % Return the final value
    left = sum;
end

function right = getRightHT(fun, s, tol)
    alpha = s +0.1;
    beta = s+25;
    
    iterates = 1; maxIterate = 400;
    sum = 0; prev = quadnc6(fun, alpha, beta, tol);
    while (iterates< maxIterate), iterates = iterates+1;
        sum = sum + prev;
        
        pBeta = beta;
        beta = beta + 2*(beta-alpha);
        alpha = pBeta;
        
        next = quadnc6(fun, alpha, beta, tol);
        if ((next/prev < tol) || (abs(prev-next)<tol)), break; end;
        prev = next;
    end;
    
    % Something gone wrong
    if (iterates >= maxIterate), warning('HILTRANSgetRightHT:MaxIterReached', 'maximum number of iterations reached - the integral seems to be singular'); end;
    
    % Return the final value
    right = sum;
end

function middle = getMiddleHT(~,~,~)
    middle = 0;
end