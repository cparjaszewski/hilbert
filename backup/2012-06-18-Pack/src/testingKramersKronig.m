function testingKramersKronig
    clear all;
    x = linspace(-5,5,1000);
    [rIn,iIn] = meshgrid(x,x);
    zIn = rIn + 1i.*iIn;
    zOut = test003(zIn);
    rOut = real(zOut);
    iOut = imag(zOut);
    rKKOut = KKreal(iOut);
    iKKOut = KKimag(rOut);
    %mesh(rIn,iIn,rOut,iOut);
    mesh(rIn, iIn, rKKOut, iKKOut);
end

function zOut = test001(zIn)
    a = real(zIn);
    b = imag(zIn);
    reOut = (a+b) ./ (b.^2 + 1);
    imOut = (a-b) ./ (a.^2 + 1);
    zOut = complex(reOut, imOut);
end


function zOut = test002(zIn)
    a = real(zIn);
    b = imag(zIn);
    reOut = (a+b) ./ (a.^2 + b.^2 + 1);
    imOut = (a-b) ./ (a.^2 + b.^2 + 1);
    zOut = complex(reOut, imOut);
end


function zOut = test003(zIn)
    a = real(zIn);
    b = imag(zIn);
    reOut = (a.^3 + 1)./(a.^5+1);
    imOut = (b.^2 + 1)./(b.^6+1);
    zOut = complex(reOut, imOut);
end

function rOut = KKreal(iIn)
    rOut = iIn;
end

function iOut = KKimag(rIn)
    iOut = rIn;
end