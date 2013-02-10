function HY = frthilbert(Y)
    %% Description:
    % This is a standard and efficient algorithm ferforming the discrete
    % Hilbert transform based on two discrete Hartley transforms in O(n log
    % n) with zero padding
   
    %% INPUT:
    % Y - an discrete array of values (ordinates)
    
    %% OUTPUT:
    % HY   - an array of Hilbert Transform values * (ordinates)
    
    %% Based on:
    % 1. Soo-Chang Pei and Sy-Been Jaw - "Computation of Discrete Hilbert
    % Transform through Fast Hartley Transform"
    % 2. http://en.wikipedia.org/wiki/Discrete_Hartley_transform
    
    %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform in~nonlinear optics"
    % krzysztof.parjaszewski@gmail.com
    
    %% The algorithm:
        
    % We perform the zero padding to the next power-of-two length
    N = max(size(Y));
    M = 2 ^ ceil(log2(N));
    Y = [Y,zeros(1,M-N)];
    
    % Discrete Hartley transform boosted-up to O(n log n)
    HF = frt(Y);
    
    % Defining the HH vector
    O1 = ones(1,floor(M/2)-1);O2 = -ones(1,ceil(M/2)-1);HH = [0, O1, 0, O2];
    
    % Defining the time reversal of HF
    TRHF = HF([1, M:-1:2]);
    
    % Based on the convolution theorem we get the Hartley-Hilbert transform
    % of X, so in the last step we need to perform the inverse Hartley
    % transform
    IHX = TRHF .* HH;
    
    % Inverse Hartley transform  boosted-up to O(n log n)
    HY = - 1/M .* frt(IHX);
    
    % The final output vector
    HY = HY(1:N);
end

