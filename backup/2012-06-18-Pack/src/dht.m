function YH = dht(X)
    %% Description:
    % HilbertTransform - Discrete Hartley transform, based on Soo-Chang Pei 
  
    %% INPUT:
    % X - an discrete vector of values (abscissas)
    
    %% OUTPUT:
    % YH   - an array of Hilbert Transform values  (ordinates)
    
    %% Based on:
    % Soo-Chang Pei  - "The Hilbert transform"
    % 2. http://en.wikipedia.org/wiki/Hermite_polynomial
    % 3. http://en.wikipedia.org/wiki/Hilbert_transform

    %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw, Summer 2011-2012]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform in~nonlinear optics"
    % krzysztof.parjaszewski@gmail.com

    N = max(size(X));
    X = reshape(X,1,N);
    M = 2 ^ (ceil(log2(N)));
    X = [X,zeros(1,M-N)];
    B = bitreversecopy(X,M);
    
    cas=@(x) (cos(x)+sin(x));
    
    for s=1:log2(M)
        m = 2^s;
        c = 1;
        for j=0:(m/2-1) 
            t = c .* B((j:m:M-1)+1+m/2);
            u = B((j:m:M-1)+1);
            B((j:m:M-1)+1) = u + t;
            B((j:m:M-1)+1 + m/2) = u - t;
            c = cas(2*pi*s/M);
        end
    end
    YH = B(1:N);
end


 function B = bitreversecopy(X,M)
    BIN = dec2bin(0:M-1);
    REV = flipdim(BIN,2);
    NUMS = bin2dec(REV);
    B = X(NUMS+1);
 end