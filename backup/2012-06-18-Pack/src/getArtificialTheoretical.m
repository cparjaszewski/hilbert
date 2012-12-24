function getArtificialTheoretical  
%% We will need:
% n,k,nm,kApproxFun (or nApproxFun or both)
%%

% function mainBody
    clc;
    display('The theoretical fit for easy complex function');
%     fun1 = @(w) (1./(300-w-5i));
%     nm = linspace(100,600,600);
%     
%     n = real(fun1(nm));
%     k = imag(fun1(nm));
    tab = GaSb;
    nm = tab(:,1);n = tab(:,2);k=tab(:,3);
    kApproxFun = @(x)GaSbLorentzianFun(x);
    %kApproxFun =  GaSbLorentzianFun;
    
    [sizexk sizeyk] = size(nm);
    l = max(sizexk,sizeyk);
    
    % CPV Method
    nCpvKK = ones(l,1);
    %kCpvKK = ones(l,1);
    rangeMax = max(nm);
    for i=1:l, x1 = nm(i); 
           nCpvFun = @(x)(x.*kApproxFun(x)./(x+x1)); nCpvKK(i) = 2./pi.*cpv(0,50*rangeMax,x1,2000,nCpvFun); 
           %kCpvFun = @(x)(real(fun1(x))./(x+x1)); kCpvKK(i) = -x1.*2./pi.*cpv(0,50*rangeMax,x1,2000,kCpvFun); 
    end
    nCpvKK = nCpvKK';
    %kCpvKK = kCpvKK';
    
    % Build-in Hilbert Transfortm
    nmH = linspace(nm(1),nm(l),1000);
    kH = kApproxFun(nmH);
    % hilbertK = imag(hilbert(n,l));
    hilbertN = 1./pi.*imag(hilbert(kH,1000));
    
    % Clenshaw-Curtis-33 method
    %kCcGK = zeros(1,l);
    nCcGK = zeros(1,l);
    
    % Adaptation of the quadgk() method
    %kQuadGK = zeros(1,l);
    nQuadGK = zeros(1,l);
    for  i=1:l, 
        t = nm(i);
        errDef = 0.005*mean(k);
        % f1 = @(x)(real(fun1(x))./(x+t));
        f2 = @(x)(x.*kApproxFun(x)./(x+t));
        
        %kQuadGK(i) = -2.*t./pi.*(quadgk(@(x)f1(x)./(x-t),0,t-errDef,'abstol',1e-14,'MaxIntervalCount',500000) + ...
        %quadgk(@(x)f1(x)./(x-t),t+errDef, inf,'abstol',1e-14,'MaxIntervalCount',500000)); % + ...
        %quadgk(@(x)(f1(t+x)-f1(t-x))./x,0,errDef,'abstol',1e-14,'MaxIntervalCount',500000);
    
        nQuadGK(i) = 2./pi.*(quadgk(@(x)f2(x)./(x-t),0,t-errDef,'abstol',1e-14,'MaxIntervalCount',500000) + ...
        quadgk(@(x)f2(x)./(x-t),t+errDef, inf,'abstol',1e-14,'MaxIntervalCount',500000)); % + ...
        %quadgk(@(x)(f2(t+x)-f2(t-x))./x,0,errDef,'abstol',1e-14,'MaxIntervalCount',500000);        
        
        % kCcGK(i) = -2.*t./pi.*cauchy(f1,t,0,inf,1e-14,errDef); 
        nCcGK(i) = 2./pi.*cauchy(f2,t,0,inf,1e-14,errDef); 
    end
    
    % final plot - N
    clf,p=plot(nm,k, ...
             nm,n, ...
             nm,3-nCpvKK, ...
             nm,3-nQuadGK, ...
             nmH,3+hilbertN, ...
             nm,3-nCcGK ... 
             );set(p,'LineWidth',3);title('n - GaSb KK'),legend('k','n','n - cpv','n - quadgk','n - hilbert','n - cc33');
         
    % final plot - K
%     figure,clf,plot(nm,k, ...
%              nm,n, ...
%              nm,kCpvKK,...
%              nm,kQuadGK, ...
%              nm,hilbertK, ...
%              nm,kCcGK ...
%              ),title('k'),legend('k','n','k - cpv','k - quadgk','k - hilbert','k - cc33');
end
    
function w = cauchy(f,t,a,b,tol,e) % t - osobliwoœæ, a,b - zasiêg, tol - tolerancja, e - jak blisko podchodzimy punktu
    f = fcnchk(f);
    w = quadgk(@(x)f(x)./(x-t),a,t-e,'abstol',tol,'MaxIntervalCount',50000) + ...
        quadgk(@(x)f(x)./(x-t),t+e, b,'abstol',tol,'MaxIntervalCount',50000);% + ...
        % quadgk(@(x)(f(t+x)-f(t-x))./x,0,e,'abstol',tol,'MaxIntervalCount',50000);
end

