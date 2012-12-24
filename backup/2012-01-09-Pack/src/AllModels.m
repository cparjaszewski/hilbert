classdef AllModels < handle
  %% Description:
    % AllModels - usual linear and nonlinear optical susceptibility models, found in common literature [1], [2]
    % shortly described in chapters 2.1-2.6
    
    %% Literature:
    % [1] Mark G. Kuzyk, Benjamin R. Anderson, Nathan J. Dawson, Sheng-Ting Hung, Wei Lu, Shiva Ramini, Jennifer L. Schei, Shoresh Shafei, Julian H. Smith, Afsoon Soudi, Szymon Steplewski, Xianjun Ye. Nonlinear optics. http://www.nlosource.com/LectureNotesBook.pdf, August 2010. Course book.
    % [2] Honam Yum and Selim M Shahriar. Pump-probe model for the kramers-kronig relations in a laser. Comments: 10 pages, 5 figures, Mar 2010.
    % [3] Robert W. Boyd. Nonlinear Optics, Third Edition. Academic Press, 3 edition, April 2008.
    % [4] Jae-Kuk Seo. Study of time-resolved measurement of intensity-dependent refractive index change in GaInAsP waveguides. Doctoral dissertation, Tokyo Institute of Technology, Department of Electrival and Electronic Engineering, 2-12-1 Ookoyama, Meguro-Ku, Tokyo, 152-8552, JAPAN, February 2006. Supervisor: Professor Tetsuya Mizumoto, Report number: 6444.
    % [5] David Crichton Hutchings, M Sheik-Bahae, D J Hagan, and E W Van Stryland. Kramers-kronig relations in nonlinear optics. Optical and Quantum Electronics, 24(1):1–30, 1992.
    
    %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw, 2011-2012]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform used to 
    % better understand and solve the Kramers-Kronig relations in nonlinear optics"
    % krzysztof.parjaszewski@gmail.com
    
    
    properties 
        all_models; % All usual models 
    end
    
     
    methods % Constructor 
        function obj = AllModels()
          obj.all_models = struct( ...
            'name' ...
            , { ...
                 'Linear model' ...
               , 'Linear model shifted' ...
               , 'Pump probe model' ...
               , 'Frequency mixing model' ...
               , 'Linear quantum pertburbative model' ...
               , 'Second order quantum perturbative model' ...
              } ...
            , 'f' ... 
            , {
                 @obj.simpllin ...
               , @obj.simpllin_shifted ...
               , @obj.simplpnp ...
               , @obj.simplmix ...
               , @obj.quantlin ...
               , @obj.quantqdr ...
               } ...
            , 'x' ...
            , {
                linspace(-40, 40, 100) ...
              , linspace(-5,   5, 100) ...
              , linspace(-18, 18, 100) ...
              , linspace(-10, 10, 100) ...
              , linspace(-12, 12, 100) ...
              , linspace(-28, 28, 1000) ...
              } ...
            , 'orientation' ... 
            , {
                 1 ...
               , 1 ...
               , 2 ...
               , 1 ...
               , 2 ...
               , 2 ...
               } ...  
            );
        end
    end
    
    methods (Access=public, Static) % Simple linear models
        function res = simpllin(x) % linear model
            res = -20./(x.*1i+1-20i) ./ (x.*1i +1+20i);
        end
        
        function res = simpllin_shifted(x) % simple linear model shifted to have center on zero
            res = -20./((x+20).*1i+1-20i) ./ ((x+20).*1i+1+20i);
        end
    end
    
    methods (Access=public, Static) % Simple nonlinear models
        function res = simplpnp(x1, x2, conf) % Complex two-argument model for pump-probe process
            %% Model configuration:
            if nargin < 2, x2 = 1.3; end
            if nargin < 3 
              G = 1;
              n0 = 1;
              eta = 1;
              theta = 1.4;
              Omega1 = 4.3;
              gammaBA = -0.1;
            else % Read configuration
              G = conf.G; 
              n0 = conf.n0; 
              eta = conf.eta;
              theta = conf.theta;
              Omega1 = conf.Omega1;
              gammaBA = conf.gammaBA; 
            end
            
            %% The value
            res = - ( (G .* n0 .* gammaBA .* 1i) ./ (x2 + x1 + 1i .* eta) .* ...
                  ( 1 - (Omega1 .^2 .* (x2 - x1 + 1i .* eta ) .* (x1 + 2 .* 1i .* eta )) ./ ...
                  ( (x2 - 1i .* eta) .* ((x1 + 1i .* theta) .* (x2 + x1 + 1i .* eta) .* (x1 - x2 + 1i .* eta ) - Omega1.^2 .* (x1 + 1i .* eta )) .* 2 ) ) );
        end

        function res = simplmix(x1, x2, conf) % Complex two-argument model for frequency mixing process
          %% Model configuration:
          if nargin<2, x2 = 1.3; end;
          if nargin<3,
              h=1;
              N=1;
              T1=1;
              T2=2;
              eps0=1;
              muBA=1;
              omega0=1;
              Omega2=1;
          else % Read configuration
              h = conf.h; 
              N = conf.N; 
              T1 = conf.T1;
              T2 = conf.T2;
              eps0 = conf.eps0;
              muBA = conf.muBA;
              omega0 = conf.omega0;
              Omega2 = conf.Omega2; 
          end
          
          %% The value is taken from [3] and described in Thesis in chapter 2.3:
          % (2.3.6)
          D = (x1 + 1 / T1) .* (x1 - x2 + 1i/T2) .* (x1 + x2 + 1i/T2) - Omega2.^2 .* (x1 + 1i/T2); 
          % (2.3.5a - description 5)
          Dc = real(D) - 1i * imag(D); 
          % (2.3.5a - numerator)
          res = (2 .* N .* omega0 .* muBA .^ 4 ) ./ (3 .* eps0 .* h .^3) .* ( -x1 - x2 -1i/T2) .* (x1 + 2i/T2); 
          % (2.3.5a - denominator)
          res = res ./ ((x2 +1i/T2) .* (x2 + x1 + 1i/T2) .* Dc); 
        end  
    end
    
    methods (Access=public, Static) % Quantum-perturbative models
        function res = quantlin(x1, conf) % Linear quantum-perturbative model    
          %% Model configuration:
          if nargin < 2, 
              conf.N  = 8;
              conf.e0 = 1.4;
              conf.h  = -2.7;
              conf.M  = [[3,-0.5];[1.2,2.4]];
              conf.O  = [-3,13];
              conf.G  = [0.7,2.3];
          end
          
          % Read configuration
          N  = conf.N;
          e0 = conf.e0;
          h  = conf.h;
          M  = conf.M;
          O  = conf.O;
          G  = conf.G; 

          %% The value:
          res = (N ./ e0 ./ h) .* ( ...
                   M(1,1) .* M(2,1) ./ (O(1) - x1 - 1i .* G(1)) ...
                 + M(2,1) .* M(1,1) ./ (O(1) + x1 + 1i .* G(1)) ...
                 + M(1,2) .* M(2,2) ./ (O(2) - x1 - 1i .* G(2)) ...
                 + M(2,2) .* M(1,2) ./ (O(2) + x1 + 1i .* G(2)) ...
              );
        end
        
        function res = quantqdr(x1, x2, conf)
            %% Model configuration:
            if nargin < 2, x2 = 1.3; end;
            if nargin < 3
              N = 5;
              e0 = 1;
              h = -1;
              M = [[1,2];[-1,-2]];
              Omega = [[3, 16];[4, 12]];
              gamma = [[1,2];[-1,3]];
            else % Read configuration
              N     = conf.N;     % for example:  5;
              e0    = conf.e0;    % for example:  1;
              h     = conf.h;     % for example: -1;
              M     = conf.M;     % for example: [[1,2];[-1,-2]];
              Omega = conf.Omega; % for example: [[3, 16];[4, 12]];
              gamma = conf.gamma; % for example: [[1,2];[-1,3]];
            end
            
            
            %% The value:
            s = 0;
            for l=1:2,
              for m=1:2,
                  for n=1:2,
                    s = s + M(l,n).*M(n,m).*M(m,l)./(Omega(n,l)-x1-x2-1i.*gamma(n,l))./(Omega(m,l)-x1-1i.*gamma(m,l)) ...
                          + M(l,n).*M(n,m).*M(m,l)./(Omega(n,l)-x1-x2-1i.*gamma(n,l))./(Omega(m,l)-x2-1i.*gamma(m,l)) ...
                          + M(l,n).*M(n,m).*M(m,l)./(Omega(m,n)-x1-x2-1i.*gamma(m,n))./(Omega(n,l)+x2+1i.*gamma(n,l)) ...
                          + M(l,n).*M(n,m).*M(m,l)./(Omega(m,n)-x1-x2-1i.*gamma(m,n))./(Omega(n,l)+x2+1i.*gamma(n,l)) ... 
                          + M(l,n).*M(n,m).*M(m,l)./(Omega(n,m)+x1+x2+1i.*gamma(n,m))./(Omega(m,l)-x1-1i.*gamma(m,l)) ...
                          + M(l,n).*M(n,m).*M(m,l)./(Omega(n,m)+x1+x2+1i.*gamma(n,m))./(Omega(m,l)-x1-1i.*gamma(m,l)) ...
                          + M(l,n).*M(n,m).*M(m,l)./(Omega(m,l)+x1+x2+1i.*gamma(m,l))./(Omega(n,l)+x1+1i.*gamma(n,l)) ...
                          + M(l,n).*M(n,m).*M(m,l)./(Omega(m,l)+x1+x2+1i.*gamma(m,l))./(Omega(n,l)+x2+1i.*gamma(n,l));
                 end;
              end;
            end;
            res = s .* N./2./e0./(h.^2);
        end
    end
end