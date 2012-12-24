function res=pertdata(omega, conf)
    %% Description:
    % Source file used for modeling the simple, linear quantum-perturbative model

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
        testhncX1; 
        res=toc;
        msg = ['Linear quantum-perturbative model test took: ' num2str(res) ' seconds'];
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

%% %% TEST ROUTINES
%% % A) htran:

function testhtran
    %% 2-dimensional test for htran method for the pertdata model
    sus1 = pertdata(x);
    plot(x,real(sus1),x,imag(sus1));
    [~,hsi] = htran(x,real(sus1));
    [~,hsr] = htran(x,imag(sus1));
    plot(x,real(sus1),x,hsr),figure,plot(x,imag(sus1),x,hsi);
    plot(x,real(sus1),x,-hsr),figure,plot(x,imag(sus1),x,hsi);
    plot(x,real(sus1),x,hsr),figure,plot(x,imag(sus1),x,-hsi);
    plot(x,real(sus1)-hsr),figure,plot(x,imag(sus1)+hsi);
    plot(x,real(sus1),x,hsr);hold on,plot(x,real(sus1)-hsr,'linewidth',4,'color','yellow');hold off;
    plot(x,real(sus1),x,hsr);hold on,plot(x,hsr-real(sus1),'linewidth',4,'color','yellow');hold off;
    figure,plot(x,imag(sus1),x,hsi);hold on,plot(x,-hsi-imag(sus1),'linewidth',4,'color','yellow');hold off;
    plot(x,imag(sus1),x,-hsi);hold on,plot(x,-hsi-imag(sus1),'linewidth',4,'color','yellow');hold off;
    plot(x,real(sus1),x,hsr,x,-hsi);hold on,plot(x,hsr-real(sus1),'linewidth',4,'color','yellow');hold off;
    figure,plot(x,imag(sus1),x,-hsi,x,hsr);hold on,plot(x,-hsi-imag(sus1),'linewidth',4,'color','yellow');hold off;
    plot(x,real(sus1),x,hsr,x,-hsi);hold on,plot(x,hsr-real(sus1),'linewidth',4,'color','yellow');hold off;
    figure,plot(x,imag(sus1),x,-hsi,x,hsr);hold on,plot(x,-hsi-imag(sus1),'linewidth',4,'color','yellow');hold off;
    figure,plot(x,imag(sus1),x,-hsi);
    plot(x,imag(sus1));
end

function testhncX1
    %% Description:
    % One-dimensional tests for linear perturbative model (fixed Delta)

    %% Quadrature configuration:
    pts = 200;
    a = -18;
    b = 18;
    n = 16;
    tol = 10^(-1);
    cs = 0.003;
    wrn = false; % we don't want the warnings here
    
    %% Model configuration:
    conf.N  = 8;
    conf.e0 = 1.4;
    conf.h  = -2.7;
    conf.M  = [[3,-0.5];[1.2,2.4]];
    conf.O  = [-3,13];
    conf.G  = [0.7,2.3];
    
    %% Main calculations:
    X = linspace(a, b, pts); % input vector
    Z = pertdata(X, conf);
    [F,H]=hncX(@(x) real(pertdata(x, conf)),a,b,tol, n, cs,pts,wrn);   
    figure,plot(X,real(Z),F,-H,F,imag(Z)),title('X,real(Z),F,-H,F,imag(Z),F,-H+imag(Z)');
    hold on, plot(F,-H-imag(Z),'Color','Yellow','LineWidth',4);
    [F1, HT1] = hncX(@(x) imag(pertdata(x, conf)),a,b,tol, n, cs,pts,wrn); 
    hold off;
    figure,plot(X,imag(Z),F1,HT1,F1,real(Z)),title('X,imag(Z),F1,HT1,F1,real(Z),F1,HT1-real(Z)');
    hold on, plot(F1,HT1-real(Z),'Color','Yellow','LineWidth',4);    
end


