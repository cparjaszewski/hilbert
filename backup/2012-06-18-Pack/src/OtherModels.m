classdef OtherModels < handle
    %% Description:
    % OtherModels - unusual nonlinear optical susceptibility models, found in literature [1], [2]
    % shortly described in chapters 2.7 and 2.8
    
    %% Literature:
    % [1] D C Hutchings, M Sheik-Bahae, D J Hagan, and E W Van Stryland. Kramers-kronig relations in nonlinear optics. Optical and Quantum Electronics, 24(1):1–30, 1992.
    % [2] Jae-Kuk Seo. Study of time-resolved measurement of intensity-dependent refractive index change in GaInAsP waveguides, PhD Thesis, Tokyo Institute of Technology,
    %     Department of Electrival and Electronic Engineering, Report number: 6444
    
    %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw, 2011-2012]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform in~nonlinear optics"
    % krzysztof.parjaszewski@gmail.com
   
    properties 
        h_models; % Models from [1]
        s_models; % Models from [2]   
    end
    
    methods % Constructor
      function obj = OtherModels()
        obj.h_models = struct('name', ...
                             {'Two Photon Absorption [1]', ...
                               'Raman [1]', ...
                               'Linear Stark [1]', ...
                               'Quadratic Stark [1]'}, ...
                            'f', ...
                              {@obj.ftpa_h, ...
                                @obj.fr_h, ...
                                @obj.fls_h, ...
                                @obj.fqs_h}, ...
                            'g', ...
                              {@obj.gtpa_h, ...
                                @obj.gr_h, ...
                                @obj.gls_h, ...
                                @obj.gqs_h} ...
                            );
         obj.s_models = struct('name', ...
                             {'Two Photon Absorption [2]', ...
                                'Raman [2]', ...
                                'Stark [2]'}, ...
                            'f', ...
                             {@obj.ftpa_s, ...
                                @obj.fr_s, ...
                                @obj.fs_s}, ...
                            'g', ...
                            {@obj.gtpa_s, ...
                              @obj.gr_s, ...
                              @obj.gs_s} ...
                            );    
        end
    end
    
    methods (Access=public, Static) % [1] Hutchings f
        % Two Photon Absorption function F_TPA(x1, x2)
        function res = ftpa_h(x1, x2)
            tpacondition = ones(size(x1));
            tpacondition((x1 + x2) < 1) = 0;
            val = (x1 + x2) .^ 3 .* (x1 + x2 - 1) .^ (3/2) .* x1.^(-3) .* x2.^(-4);
            res = Utils.getres(tpacondition, val);
        end

        % Raman function F_R(x1, x2)
        function res = fr_h(x1, x2)
            rcondition = ones(size(x1));
            rcondition((x1 - x2) < 1) = 0;
            val = (x1 - x2) .^ 3 .* (x1 - x2 - 1) .^ (3/2) .* x1.^(-3) .* x2.^(-4);
            res = Utils.getres(rcondition, val);
        end

        % Linear Stark function F_LS(x1, x2)
        function res = fls_h(x1, x2)
            scondition = ones(size(x1));
            scondition(max(x1, x2) < 1) = 0;
            val = -(x1-1) .^ (3/2) ./ (64 .* x1 .* x2 .^ 4);
            res = Utils.getres(scondition, val);
        end

        % Quadratic Stark function F_QS(x1, x2)
        function res = fqs_h(x1, x2) 
            qcondition = ones(size(x1));
            qcondition(max(x1, x2) < 1) = 0;
            val = - (1 ./ (x1 - x2) + 1 ./ (x1 + x2)) ...
                 ./ ( 1024 .* x1 .* x2 .^ 2 .* (x1 - 1) .^ (1/2)); 
            res = Utils.getres(qcondition, val);
        end
    end
    
    methods (Access=public, Static) % [1] Hutchings g
      % Two Photon Absorption function G_R(x1, x2)
      function res = gtpa_h(x1, x2)
            tpacondition = ones(size(x1));
            tpacondition((x1 + x2) < 1) = 0;
            val = (x1 + x2) .^ 3 .* (x1 + x2 - 1) .^ (3/2) .* x1.^(-3) .* x2.^(-4);
            res = Utils.getres(tpacondition, val);
        end

        % Raman function G_R(x1, x2)
        function res = gr_h(x1, x2)
            rcondition = ones(size(x1));
            rcondition((x1 - x2) < 1) = 0;
            val = (x1 - x2) .^ 3 .* (x1 - x2 - 1) .^ (3/2) .* x1.^(-3) .* x2.^(-4);
            res = Utils.getres(rcondition, val);
        end

        % Linear Stark function G_LS(x1, x2)
        function res = gls_h(x1, x2)
            scondition = ones(size(x1));
            scondition(max(x1, x2) < 1) = 0;
            val = -(x1-1) .^ (3/2) ./ (64 .* x1 .* x2 .^ 4);
            res = Utils.getres(scondition, val);
        end

        % Quadratic Stark function G_LS(x1, x2)
        function res = gqs_h(x1, x2) 
            qcondition = ones(size(x1));
            qcondition(max(x1, x2) < 1) = 0;
            val = -(1 ./ (x1 - x2) + 1 ./ (x1 + x2)) ./ ( 1024 .* x1 .* x2 .^ 2 .* (x1 - 1) .^ (1/2)); 
            res = Utils.getres(qcondition, val);
        end
      
    end
    
    methods (Access=public, Static) % [2] SEO f
        function res = ftpa_s(x1, x2)
          tpacondition = ones(size(x1));
          tpacondition((x1 + x2) < 1) = 0;
          val = (x1 + x2) .^ 3 .* (x1 + x2 - 1) .^ (3/2) .* x1 .^ (-3) .* x2 .^ (-4);
          res = Utils.getres(tpacondition, val);
        end
        
        function res = fr_s(x1, x2)
          rcondition = ones(size(x1));
          rcondition((x1 - x2) < 1) = 0;
          val = (x1 - x2) .^ 3 .* (x1 - x2 - 1) .^ (3/2) .* x1 .^ (-3) .* x2 .^ (-4);
          res = Utils.getres(rcondition, val);
        end
        
        function res = fs_s(x1, x2)
          scondition = ones(size(x1));
          scondition(x1 < 1) = 0;
          val = - 4 .* (x1 - 1) .^ (3/2) ./ (x1 .^2 .* (x1 + x2) .* (x1 - x2)) ...
                - 4 .* (x1 - 1) .^ (1/2) ./ (x1 .* (x1 + x2) .* (x1 - x2)) ...
                - (x1-1) .^ (-1/2) ./ ((x1 + x2) .* (x1 - x2)) ...
                - 2 .* (x1 - 1) .^ (3/2) .* x2 .^ (-4) ...
                - 6 .* (x1 - 1) .^ (3/2) .* (x1) .^ (-2) .* (x2) .^ (-2) ...
                - 9 .* (x1 - 1) .^ (1/2) .* x1 .^ (-1) .* x2 .^ (-2) ...
                - 3/4 .* (x1 - 1) .^ (-1/2) .* x2 .^ (-2);
          res = Utils.getres(scondition, val);
        end
        
        
    end
    
    methods (Access=public, Static) % [2] SEO g
        function res = gtpa_s(x, y)
          res = t2(x, y) + t2(-x, y);
        end
        
        function res = gr_s(x, y)
          res = t2(x, -y) + t2(-x, -y);
        end
        
        function res = gs_s(x, y)
          res = s2(x, y) + s2(-x, y) + s2(x, -y) + t2(-x, -y);
        end
        
        % Auxiliary function T_2
        function res = t2(x, y)
            res = (x+y).^3 .* (1-x-y).^(3/2)./(x.^4 .* y.^4) ...
                  - (1./(x .^ 4 .* y) ... 
                  + 3./(x .^ 2 .* y .^ 3)) .* (1-y).^(3./2) ...
                  + 9./(2.*x.^2.*y.^2) .* (1-y).^(1/2) ...
                  - (3./8)./(x.^2 .* y.^2) .* (1-y).^(1/2) ...
                  - (3./8)./(x.^2.*y) .* (1-y).^(1/2);
        end
        
         % Auxiliary function S_2
        function res = s2(x, y)
            res =  - (1 ./ (x .* y .^ 4) + 3./(x .^ 3 .* y .^ 2)) .* (1 - x) .^ (3/2) ...
                   + 9 ./ ( 2 .* x .^ 2 .* y .^ 2) .* (1 - x) .^ (1/2) ...
                   - 3./(8 .* x .* y .^ 2) .* (1- x) .^ (-1/2) ...
                   - 4 ./( x .^ 2 .* y .^ 2) ...
                   + 1 ./ (x .^ 2 - y .^ 2) .* (2 ./ y .^ 3 .* ( 1 - y) .^ (3/2) ...
                   - 2 ./ x .^ 3 .* ( 1 - x).^3/2 ...
                   - 2 ./ y.^2 .* (1 - y).^(1/2) ...
                   + 2 ./ (1 - x) .^ (1/2) ...
                   + 1 ./ (2 .* y) .* (1 - y) .^ (1/2) ...
                   - 1 ./ (2 .* x) .* (1 - x) .^ 2);
        end
    end
    
end

