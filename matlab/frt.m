function res = frt(X)
    %% Description:
    % Fast Hartley Transform algorithm. This algorithm is faster than simple FFT because it works only in
    % real domain, which involves half the number of multiplications
      
    %% INPUT:
    % X - an discrete array of values (ordinates)
    
    %% OUTPUT:
    % HY   - an array of Discrete Hartley Transform values * (ordinates)
    
    %% Based on:
    % 1. Ronald F. Ullman - "An algorithm for the Fast Hartley Transform"
    % 2. Krzysztof Loryœ - polish notes to the lecture of "Algorithms and
    % Data Structures" held in the Institute of Computer Science,
    % University of Wroclaw
  
	  %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform in~nonlinear optics"
    % krzysztof.parjaszewski@gmail.com
    
    %% The algorithm:
    N = length(X);
    
    % We precalculate the CAS table for vectors of length less or equal 8:
    CAS = cell(8);
    CAS{1} = cas(0);
    for k=2:8, CAS{k} = cas(2*pi*((0:(k-1))'*(0:(k-1)))/k); end
    
    % Now we are sure, that X vector is of size M which is a power value of 2.
    DHT = dofrt(X,N,CAS);
    
    % The final vector must be truncated because it was arbitrary enlarged
    % to be with length of power of two
    res = DHT(1:N);
end

function HY = dofrt(X, N, CAS)
    %% Recursive, divide & conquer radix-2 Fast Hartley Transform
    
    %% INPUT:
    % X   - an discrete array of values (ordinates)
    % D   - the length of X vector
    % CAS - the precalculated cell of 1:8-8 matrix for performing the instant FRT
    
    %% OUTPUT:
    % HY   - an array of Discrete Hartley Transform values * (ordinates)
    
    % Fast Hartley Trasform for N lower than 8
    if N<=8, 
        HY=X * CAS{N} ; % Instant FRT for 8-element array
        return;
    end;
    % We split the input vector into two equal vectors
    X1 = X(1:2:end); X2 = X(2:2:end); % In MATLAB first index equals 1, not 0 - but in literature X1 is called even, X2 - odd indices.
    
    % Divide & conquer approach
    HT1 = dofrt(X1,N/2,CAS); %FRT of even
    HT2 = dofrt(X2,N/2,CAS); %FRT of odd
    
    % We precalculate the lower and upper part of the output array 
    low    = 1:1:N/2;
    revlow = [1, N/2:-1:2];
    arg    = (2*pi/N) .* (low-1);
    CS     = cos(arg) .* HT2(low) + sin(arg) .* HT2(revlow);
    H      = HT1(low);
    
    % The lower and upper part of DHT is combined with the following
    % radix-2 way
    HY = [H + CS, H - CS];
end

function res = cas(X), res = cos(X) + sin(X); end