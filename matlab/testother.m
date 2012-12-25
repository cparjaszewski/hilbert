function testother
    testhilbert; return;
    x = linspace(-10,10,600);
    [x1,x2] = meshgrid(x);
    
    %FTPA:
    res = ftpa(x1, x2);
    
    %GTPA:
    gres = t2(x1,x2) + t2(-x1,x2);
    
    % HILBERT TRANSFORM:
    hvals = zeros(length(x));  iter = 0; l = length(x); h = waitbar(0,'Please wait...');
    for i = x
        iter = iter + 1;
        fivals = res(iter,:);
        hvals(iter, :) = hfthilbert(fivals);
        waitbar(iter/l);
    end, close(h);
    
    %CALC
    logres = dolog(res);    logres(isnan(logres))   = 0;
    logigres = dolog(imag(gres));
    logrgres = dolog(real(gres));
     loggres = dolog(gres);  loggres(isnan(loggres))   = 0;
    
    %PLOT
    mesh(x1,x2,logres)         , title('log of the res');
%     figure,mesh(x1,x2,res)     , title('res');  
    figure,mesh(x1,x2,logigres) , title('log of the imag of the gres');
    figure,mesh(x1,x2,logrgres) , title('log of the real of the gres');
    figure,mesh(x1,x2,imag(gres))    , title('imag gres');
    figure,mesh(x1,x2,real(gres))    , title('real gres');
    figure,loghvals = dolog(hvals);mesh(x1,x2,loghvals),title('log of the hvals');
    
end

function output = dolog(input)
     output = zeros(size(input));
     upper = input; upper(input < 0) = 0;
     lower = input; lower(input > 0) = 0; lower = - lower;
     output  = output + log(upper+1);
     output  = output - log(lower+1); 
end

function res = t2(x, y)
    res = (x+y).^3 .* (1-x-y).^(3/2)./(x.^4 .* y.^4) - (1./(x .^ 4 .* y) ... 
        + 3./(x .^ 2 .* y .^ 3)) .* (1-y).^(3./2) + 9./(2.*x.^2.*y.^2) .* (1-y).^(1/2) ...
        - (3./8)./(x.^2 .* y.^2) .* (1-y).^(1/2) ...
        - (3./8)./(x.^2.*y) .* (1-y).^(1/2);
end

function res = ftpa(x1, x2)
    tpacondition = ones(size(x1));
    tpacondition((x1 + x2) < 1) = 0;
    val = (x1 + x2) .^ 3 .* (x1 + x2 - 1) .^ (3/2) .* x1 .^ (-3) .* x2.^(-4);
    res = tpacondition .* val;
    res(isnan(res) | (isfinite(res)==0)) = 0;
end

function testhilbert
     len = 440;
     x = linspace(-10,10,len);
     [X1, X2] = meshgrid(x);
     funZ2 = @ftpa;
     Z = ftpa(X1,X2);
     inh = waitbar(0,'Please wait','Name','Hilbert Transform with Clenshaw-Curtis');
     H = zeros(size(X1));h1 = waitbar(0,'Please Wait...');sumtime = 0;
     for i = 1:1:len, 
       tStart=tic;
       hfunZ2 = @(x) funZ2(x,X2(i)); H(i,:) = htrancc(hfunZ2, x, 10^(-3), inh)'; 
       tElapsed=toc(tStart);sumtime = sumtime + tElapsed;
       timeleft = sumtime*(len-i)/i; minleft = floor(timeleft/60); secleft = timeleft - minleft*60;
       waitbar(i/len,h1,['General progress: ', num2str(i/len*100),'%, time left: ', num2str(minleft), 'm ',num2str(secleft),'s' ]);
    end; close(h1);close(inh);
    
    figure,mesh(X1,X2,dolog(H));
    figure,mesh(X1,X2,dolog(Z));
end