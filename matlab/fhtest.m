function fhtest(fun,a,b,m,funh,p)
    if nargin<6,p=1000;end;
    XF = linspace(a,b,p);
    YF = fun(XF);
    doccN = 1024;
    % We change the integration interval from (a,b) to (-1,1) to make it
    % proper for the Clenshaw-Curtis Quadrature rutine:
    const = pi;
    fun11 = @(t) (const * fun((t./2).*(b-a)+b+a));
    funcos11 = @(p,t) (const * fun((t./2).*(b-a)+b+a).*cos(p.*((t./2).*(b-a)+b+a)));
    funsin11 = @(p,t) (const * fun((t./2).*(b-a)+b+a).*sin(p.*((t./2).*(b-a)+b+a)));
    fsum = 0;
    for k = 0:1:m
        if k<1, an = 1/(2*pi) * docc(fun11,doccN ); else an = 1/pi*docc(@(t)funcos11(k,t),doccN ); end;
        bn = 1/pi * docc(@(t)funsin11(k,t),doccN ); 
        fsum = fsum + an .* cos(k*XF) + bn .* sin(k*XF);
    end
    
    hsum = 0;
    for k = 1:1:m
        han = 1/pi*docc(@(t)funsin11(k,t),doccN ); 
        hbn = 1/pi*docc(@(t)funcos11(k,t),doccN ); 
        hsum = hsum + han .* cos(k*XF) + hbn .* sin(k*XF);
    end
    if nargin< 5, plot(XF,YF,XF,fsum,XF,hsum); else plot(XF,YF,XF,fsum,XF,hsum,XF,funh(XF)); end;
end