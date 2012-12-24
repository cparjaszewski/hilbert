function val = docc(fun, N, a, b)
    %% Description:
    % This is one point algorithm and should be used in algorithms that
    % need to determine the full vector values
    % Look precisely:
    % abcde -> abcdedcb ( there is no second 'a' )
    
    %% INPUT:
    % fun - an function_handle 
    % N - number of points to operate with tolerance (scalar) [optional, default:256]
    % a - the lower integration bound (scalar) [optional, default: -1]
    % b - the upper integration bound (scalar) [optional, default:  1]
    
    %% OUTPUT:
    % val  - value of integral
    
    %% Based on:
    % 1. http://fourier.eng.hmc.edu/e161/lectures/dct/node2.html
    % (unfortunatelly this is DCT-II, so another type)
    % 2. http://en.wikipedia.org/wiki/Clenshaw–Curtis_quadrature
    % The very important quote: 
    % "For example, a DCT-I of N=5 real numbers abcde is exactly equivalent
    % to a DFT of eight real numbers abcdedcb (even symmetry), divided by
    % two."
    
    %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw, 2011-2012]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform in~nonlinear optics"
    % krzysztof.parjaszewski@gmail.com
    
    %% The algorithm:
    if nargin<4, b =   1; end;
    if nargin<3, a =  -1; end;
    if nargin<2, N = 256; end;

    % Mapping [a, b] into [-1, 1]
    C = (b+a)/2;
    T = (b-a)/2;
    
    % A new function to be integrated from -1 to 1    
    newfun = @(x) fun((x .* T ) + C);

    % Getting the value with standard -1,1 Clenshaw Curtis integration routine
    val = onetoonecc(newfun,N);
    
end
    
function val = onetoonecc(fun, N)
    %% CC integration from -1 to 1:
    
    N2  = N/2;                                           % We make this division once, for faster computation
    ps  = 0:1:N2;                                        % Array of N/2+1 points (where N2 = N/2)
    xn  = cos(ps.*pi./N);                                % N2+1 abscissaas from the 0..pi range
    yn  = fun(xn) + fun(-xn);                            % N2+1 ordinates for fun(cos(n*pi/N))
    yn(isnan(yn) | (isfinite(yn)==0)) = 0;               % We shall set all NaNs and Infs to zero
    g   = real(fft(yn(1+[0:N2 (N2-1):-1:1])))./2;        % Look at quote in the base description about the connection of DCT-I and DFT/FFT
    cn  = [g(1), g(2:N2)+g(N:-1:N2+2), g(N2+1)];         % We take symmetrical values and add them together (ommiting first and last)
    dn  = [2, 2./(1-(2.*(1:1:(N2-1))).^2), 2./(1-N^2)];  % This is the final vector for multiplying
    val = (dn * cn' )/ N;                                % Final cc integral calculation
end