function I = cc1(fun,n)                          % Kwadratura C-C - (n+1)-punktowa
    f  = fcnchk(fun);                            % Akceptuj rozne formy funkcji

    x  = cos(pi*(0:n)'/n);                       % wezly Czebyszewa
    fx = feval(f,x)/(2*n);                       % Wartosci f w tych wezlach
    g  = real(fft(fx([1:n+1 n:-1:2])));          % Szybka transformacja Fouriera
    a  = [g(1); g(2:n)+g(2*n:-1:n+2); g(n+1)];   % Wspolczynniki Czebyszewa

    w  = 0*a'; w(1:2:end) = 2./(1-(0:2:n).^2);   % Calki z wiel. Czebyszewa
    % I = w*a;                                   % Calka funkcji f
    I  = w(end:-1:1)*a(end:-1:1);                % Tak lepiej sumowac!
end