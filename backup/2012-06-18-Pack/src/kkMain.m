function kkMain()
% What is to test and compare:
% 1. Cauchy Principal Value
% 2. QUADGK()
% 3. Schemat C-C
% 4. Differential equations (set of 8 different methods EmRK, ERK)
% 5. Keller's one

% How to test and compare:
% a. different approximation schemes
% b. noise
% c. different approximation parameters

% Different Kramers Kronig Models:
% linear
% nonlinear - version 1
% nonlinear - version 2
% nonlinear - version 3

% Some other models - think-about
    
    kk = KramersKronig(@fun001);
    
    kk.runCPV;
    kk.runQuadkGK
    kk.runCC;
    kk.runRK;
    kk.runKeller;
end

function z = fun001(w)
    z = exp(-i * w);
end