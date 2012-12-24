function h = hcc(fun, y, tol)
    % We need to change variable in the infinite integration x =
    % t/(1-t.^2) and we assume limit(fun(t)/t)=0 for +-infinity as should
    % be in Hilbert transform
    ccFun = @(t) (abs(t)<1.*(-fun(t./(1-t.^2))+fun(y)).*(1+t.^2)./(1-t.^2)./(-y+y.*t.^2+t));
    
    % We set N to the initial value:
    N = 8; lastVal = 0;
    while (true)
        h =  1/pi*docc(ccFun,N);
        
        % A simple adaptive method - if the relative error is to large - we
        % double the points count:
        if (abs(lastVal-h)/(abs(h)+tol)<tol),break; end;
        N = N*2;
        lastVal = h;
    end
end