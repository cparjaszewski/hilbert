function [HY, Y] = fourhtrans(fun, X, m, doccN)
    %% Description: 
    % This is a hilbert transform function based on the Fouries series
    % approximation of the given function.
    
    %% INPUT:
    % fun - an discrete array of values        (function handle)
    % X - interval to perform Hilbert transform  (abscissas)
    % m - order of Fourier series expantions   (scalar) [optional, default:30]
    % doccN - number of points for CC integral (scalar) [optional, default:1024]
    
    %% OUTPUT:
    % HY - an array of Hilbert Transform values (ordinates)
    % Y  - an array of fun(X) approximation     (ordinates)
    
    %% Based on:
    % 1. M. Johansson - "The Hilbert transform"
    % 2. http://en.wikipedia.org/wiki/Fourier_series
    % 3. http://en.wikipedia.org/wiki/Hilbert_transform
    
    %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw, Summer 2011]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform used to 
    % better understand and solve the Kramers-Kronig relations in nonlinear optics"
    % krzysztof.parjaszewski@gmail.com
    
    %% The algorithm:
    if nargin<3,m=30; end;
    if nargin<4,doccN=1024; end;

    % First we calculate the X boundaries
    a = min(X);       
    b = max(X);
    n = length(X);

    % Then we enlarge 5 times the investigated interval:
    down = a-(b-a)*2; 
    up   = b+(b-a)*2;
    XX = linspace(down, up, 5*n);

    % We measure the important values
    T = (up-down);       % period
    L = T/2;             % halfperiod

    % Calculation of the first component
    Yp  = 1/2 .* docc(fun, doccN, down, up);
    HYp = 0;

    % Main loop:
    for k=1:1:m
        an  = docc(@(x) fun(x) .* cos(k*pi/L .* x), doccN, down, up);
        bn  = docc(@(x) fun(x) .* sin(k*pi/L .* x), doccN, down, up);

        cs = cos(k*pi/L .* XX);
        sn = sin(k*pi/L .* XX);
        
        % Function approximation (on the NX interval)
        Yp  =  Yp + an .* cs + bn .* sn;
        % Hilbert transform approximation (on the NX interval)
        HYp = HYp + bn .* cs + an .* sn;
    end

    % Finally, we extract the [a, b] part from the [up, down] interval
    INN = (2*n+1):1:3*n;
    Y  =  Yp(INN); 
    HY = HYp(INN); 

    plot(X, fun(X), X, Y, X, HY);
end