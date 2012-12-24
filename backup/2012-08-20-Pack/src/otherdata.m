function otherdata(omega, conf)
    %% Description:
    % Source file used for modeling the other nonlinear models defined in
    % Chapters 2.7 and 2.8

    %% INPUT:
    % omega - an discrete array of values  (vector)
    % conf - model configuration           (special)
    
    %% OUTPUT:
    % res - comple array of model values   (vector)

    %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw, 2011-2012]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform in~nonlinear optics"
    % krzysztof.parjaszewski@gmail.com

    %% Test case:
    % if no arguments are given, tests are run
    if (nargin==0)
        tic
%         runtest2d('htran'); 
        runtest3d('htran'); 
%         runtestdiv;
        res=toc;
        msg = ['Nonlinear other models'' test took: ' num2str(res) ' seconds'];
        display(msg);
    else
        %% The model:
        N  = conf.N;
        e0 = conf.e0;
        h  = conf.h;
        M  = conf.M;
        O  = conf.O;
        G  = conf.G; 

        %% The value
        res = (N ./ e0 ./ h) .* (...
                 M(1,1) .* M(2,1) ./ (O(1) - omega - 1i .* G(1)) ...
               + M(2,1) .* M(1,1) ./ (O(1) + omega + 1i .* G(1)) ...
               + M(1,2) .* M(2,2) ./ (O(2) - omega - 1i .* G(2)) ...
               + M(2,2) .* M(1,2) ./ (O(2) + omega + 1i .* G(2)) ...
              );
    end;
end

function runtestdiv
  x = -0.5:0.001:0.5;
 
  DVAL = divergence(x);
  
  figure;
  subplot(2, 1, 1); plot2dlog(x, imag(DVAL)); title('imag divergence')
  subplot(2, 1, 2); plot2dlog(x, real(DVAL)); title('real divergence')
end

function runtest2d(method)
  % One dimensional
  data = -0.5:0.001:0.5;
  chosenpoint = data(ceil(length(x)/1.4));
  % Set all three axes to be log scale
  set(gca, 'XScale', 'linear', 'YScale', 'linear', 'ZScale', 'linear');
  om = OtherModels();
  
  for model = om.h_models
    plotfun2d(data, model, method, chosenpoint);
  end
end

function runtest3d(method)
  % One dimensional
  data = linspace(-4.5,4.5,200);
  
  % Set all three axes to be log scale
  set(gca, 'XScale', 'linear', 'YScale', 'linear', 'ZScale', 'linear');
  om = OtherModels();
  
  for model = om.h_models
      plotfun3d(data, model, method);
  end
  
  for model = om.s_models
      plotfun3d(data, model, method);
  end
end

function plotfun2d(data, model, method, chosenpoint)
    ffun = model.f; gfun = model.g; pname = model.name; x = data;
    FVAL2D = ffun(chosenpoint, x); FVAL2D = isfinite(FVAL2D) .* FVAL2D; FVAL2D = isnumeric(FVAL2D) .* (FVAL2D);

    figure;
    subplot(2, 3, 1); plot2dlog(x,imag(FVAL2D)); title(['imag f -  ' pname ' at point: ' num2str(chosenpoint) ])
    subplot(2, 3, 2); plot2dlog(x,real(FVAL2D)); title(['real f -  ' pname ' at point: ' num2str(chosenpoint) ])
    
    switch(methodname)
        case 'htran'
            grvals = htran(imag(FVAL2D), x);
            givals = htran(real(FVAL2D), x);
        case 'hncX' 
            % even in 2d- lots of minutes of waiting
            %% Quadrature configuration:
            pts = length(x);a = min(x);b = max(x);
            n = 8; tol = 10^(2); cs = 0.03;
            wrn = false; % we don't want the warnings here
            
            %% Main calculations:
            [~,grvals] = hncX(@(y) imag(ffun(chosenpoint, y)),a,b,tol, n, cs,pts,wrn);   
            [~,givals] = hncX(@(y) real(ffun(chosenpoint, y)),a,b,tol, n, cs,pts,wrn);   

        case 'hfthilbert'            
            grvals = hfthilbert(imag(FVAL2D));
            givals = hfthilbert(real(FVAL2D));
        otherwise
           error(['Test method: ' method ' not implemented yet']); 
    end
    subplot(2, 3, 3); plot2dlog(x, real(grvals)); title(['real hilbert 2d grvals -  ' pname ' at point: ' num2str(chosenpoint) ])
    subplot(2, 3, 4); plot2dlog(x, imag(grvals)); title(['imag hilbert 2d grvals -  ' pname ' at point: ' num2str(chosenpoint) ])
    subplot(2, 3, 5); plot2dlog(x, real(givals)); title(['real hilbert 2d givals -  ' pname ' at point: ' num2str(chosenpoint) ])
    subplot(2, 3, 6); plot2dlog(x, imag(givals)); title(['imag hilbert 2d givals -  ' pname ' at point: ' num2str(chosenpoint) ])
    figure; plot(x, gfun(x)); title(['gfun -  ' pname]);
end

function plotfun3d(data, model, method)
    x = data;
    [X1, X2] = meshgrid(x, x);
    ffun = model.f; gfun = model.g; pname = model.name;
    FVAL = ffun(X1, X2); FVAL = onlynums(FVAL);
%   GVAL = gfun(X1, X2); GVAL = isfinite(GVAL) .* GVAL .* isnumeric(GVAL);
    
    if (strcmp(pname, 'Two Photon Absorption [1]') == 0), figure; end
    subplot(2, 3, 1); plot3dlog(X1, X2, imag(FVAL)); title(['imag f -  ' pname ])
    subplot(2, 3, 2); plot3dlog(X1, X2, real(FVAL)); title(['real f -  ' pname ])
%   subplot(2, 4, 3);plot3dlog(X1,X2,imag(GVAL));title(['imag g -  ' pname ])
%   subplot(2, 4, 4);plot3dlog(X1,X2,real(GVAL));title(['real g -  ' pname ])

    hrvals = zeros(length(x));  hivals = zeros(length(x));
    iter = 0;
    h = waitbar(0, 'Please wait...'); l = length(x);

    switch(method)
        case 'htran'
            for i = x
                iter = iter + 1;
                fivals = FVAL(:, iter); fivals = onlynums(fivals);
                hrvals(:, iter) = htran(imag(fivals), x)';
                hivals(:, iter) = htran(real(fivals), x)';
                waitbar(iter/l);
            end    
        case 'hncX'
            %% Extremely long (hours of hours...)
            
            %% Quadrature configuration:
            pts = length(x);a = min(x);b = max(x);
            n = 8; tol = 10^(2); cs = 0.03;
            wrn = false; % we don't want the warnings here

            for i = x
                iter = iter + 1;
                %% Main calculations:
                [~,IN] = hncX(@(y) imag(ffun(y,i)),a,b,tol, n, cs,pts,wrn);   
                hrvals(:, iter) = IN';
                [~,IN] = hncX(@(y) real(ffun(y,i)),a,b,tol, n, cs,pts,wrn);   
                hivals(:, iter) = IN';
                waitbar(iter/l);
            end 
         case 'hfthilbert'
            for i = x
                iter = iter + 1;
                fivals = ffun(x, i); fivals = onlynums(fivals);
                hrvals(:, iter) = hfthilbert(imag(fivals))';
                hivals(:, iter) = hfthilbert(real(fivals))';
                waitbar(iter/l);
            end 
        otherwise
            error(['Test method: ' method ' not implemented yet']);
    end
    
    close(h);
    subplot(2, 3, 3);plot3dlog(X1, X2, real(hrvals));title(['real hrvals 3d -  ' pname])
    subplot(2, 3, 4);plot3dlog(X1, X2, imag(hrvals));title(['imag hrvals 3d -  ' pname])
    subplot(2, 3, 5);plot3dlog(X1, X2, real(hivals));title(['real hivals 3d -  ' pname])
    subplot(2, 3, 6);plot3dlog(X1, X2, imag(hivals));title(['imag hivals 3d -  ' pname])
    %   figure;
    %   subplot(1,2,1);mesh(X1,X2,imag(FVAL) -hrvals );title(['imag 3d error -  ' pname])
    %   subplot(1,2,2);mesh(X1,X2,real(FVAL) -hivals );title(['real 3d error -  ' pname])
    % GVAL = gfun(x);
    % figure;plot2dlog(x, GVAL);title(['hilbert g -  ' pname])
end
 
%% Biargumental Models:

% Two Photon Absorption function F_TPA(x1, x2)
function res = ftpa(x1, x2)
    tpacondition = ones(size(x1));
    tpacondition((x1 + x2) < 1) = 0;
    val = (x1 + x2) .^ 3 .* (x1 + x2 - 1) .^ (3/2) .* x1.^(-3) .* x2.^(-4);
    res = getres(tpacondition, val);
end

% Raman function F_R(x1, x2)
function res = fr(x1, x2)
    rcondition = ones(size(x1));
    rcondition((x1 - x2) < 1) = 0;
    val = (x1 - x2) .^ 3 .* (x1 - x2 - 1) .^ (3/2) .* x1.^(-3) .* x2.^(-4);
    res = getres(rcondition, val);
end

% Linear Stark function F_LS(x1, x2)
function res = fls(x1, x2)
    scondition = ones(size(x1));
    scondition(max(x1, x2) < 1) = 0;
    val = -(x1-1) .^ (3/2) ./ (64 .* x1 .* x2 .^ 4);
    res = getres(scondition, val);
end

% Quadratic Stark function F_LS(x1, x2)
function res = fqs(x1, x2) 
    qcondition = ones(size(x1));
    qcondition(max(x1, x2) < 1) = 0;
    val = -(1 ./ (x1 - x2) + 1 ./ (x1 + x2)) ./ ( 1024 .* x1 .* x2 .^ 2 .* (x1 - 1) .^ (1/2)); 
    res = getres(qcondition, val);
end

% Removes nans and infinites
function res = getres(condition,val)
    res = condition .* val;
    res = onlynums(res);
end

function res = divergence(x)
    res = 1./(2*x).^6 .* (-2 - 35/8*x.^2 +1/8*x.*(3*x - 1).*(1-x).^(-1/2) - 3*x.*(1-x).^(1/2) + (1-x).^(3/2));
end

%% Hilbert transforms:

% Modified Hilbert transform (by x2) of the Two Photon Absorption function G_TPA(x1, x2)
function res = gtpa(x1, x2) 
    l1 = log(1./(x1+x2-1));
    l2 = log(-1./(x1+x2-1));
    res = 1/2.*(x1+x2-1).^(1/2).*(x2>=0).*( ...
                  l1.*x1.^4        -    l2.*x1.^4       -    l2.*x2.^4        ...
             +    l1.*x2.^4        + 4.*l1.*x1.*x2.^3   - 6.*l2.*x1.^2.*x2.^2 ...
             + 6.*l1.*x1.^2.*x2.^2 - 4.*l2.*x1.^3.*x2   +    l2.*x1.^3        ...
             + 4.*l1.*x1.^3.*x2    -    l1.*x1.^3       -    l1.*x2.^3        ...
             +    l2.*x2.^3        - 4.*l2.*x1.*x2.^3   + 3.*l2.*x1.*x2.^2    ...
             - 3.*l1.*x1.*x2.^2    + 3.*l2.*x1.^2.*x2   - 3.*l1.*x1.^2.*x2)   ...
       ./pi./x1.^3./x2.^5;
end

% Modified Hilbert transform (by x2) of the Raman function F_R(x1, x2)
function res = gr(x1, x2) 
    res = -1/128*(x1+x2-1)^(3/2)*(-(ln(-1/(x1+x2-1))*x1^2-ln(1/(x1+x2-1))*x1^2+2*ln(-1/(x1+x2-1))*x2*x1-2*ln(1/(x1+x2-1))*x2*x1+ln(-1/(x1+x2-1))*x2^2-ln(1/(x1+x2-1))*x2^2)*(-1+Heaviside(x2)))/Pi/x2^5/x1^3; 
end

% Modified Hilbert transform (by x2) of the Linear Stark function F_LS(x1, x2)
function res = gls(x1, x2) 
    res = -1/64*I*(x1-1)^(3/2)*(2*Heaviside(x2)-1)/x1/x2^5;
end

% Modified Hilbert transform (by x2) of the Quadratic Stark function F_LS(x1, x2)
function res = gqs(x1, x2) 
    res = -1/512*I*(-2*Heaviside(x2)*x1^3+2*Heaviside(x1)*x2^3+x1^3-x2^3)/x2^3/(x1-1)^(1/2)/(-x1^2+x2^2)/x1^3;
end

function plot3dlog(x,y,z)
     zvalue = dolog(z);
     % 2 is a semi-arbitrary "fudge factor", designed to make absolutely sure
     % that none of the matrix elements even come close to 0
     mesh(x, y, zvalue);
end

function plot2dlog(x,y)
    y = y .* isfinite(y);y = y .* isnumeric(y);
    plot(x,log(abs(y)).*sign(y));
%     plot(x,y);
end

function res = t2(x, y)
    res = (x+y).^3 .* (1-x-y).^(3/2)/(x.^4 .* y.^4) - (1/(x .^ 4 .* y) + 3/(x .^ 2 .* y .^ 3)) .* (1-y).^(3/2) + 9./(2.*x.^2.*y.^2) .* (1-y).^(1/2) - (3/8)./(x.^2 .* y.^2) .* (1-y).^(1/2) -(3/8)./(x.^2.*y).*(1-y).^(1/2);
end

% Makes logarithm of the positive and negative values separately
function output = dolog(input)
     output = zeros(size(input));
     upper = input; upper(input < 0) = 0;
     lower = input; lower(input > 0) = 0; lower = - lower;
     output  = output + log(upper+1);
     output  = output - log(lower+1);
     output = onlynums(output);
end

% Zeros NaN and infinities
function res = onlynums(input)
    res = input;
    res(isnan(res) | (isfinite(res)==0)) = 0;
end