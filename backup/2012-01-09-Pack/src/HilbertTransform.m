classdef HilbertTransform < handle
    %% Description:
    % HilbertTransform - model for testing all gathered Hilbert transform implementations 

    %% Used Hilbert Transforms:
    % a. htran      - based on paper by I.J. Weinbert             [MSc Thesis, chapter 4]
    % b. hncX       - based on Newton-Cotes quadrature            [MSc Thesis, chapter 5]
    % c. htrancc    - based on the Clenshaw-Curtis quadrature     [MSc Thesis, chapter 6]
    % d. hfthilbert - based on the fast Hartley transform         [MSc Thesis, chapter 7]
    % e. herhtrans  - based on Hermite transform                  [MSc Thesis, chapter 8]
    % f. fourhtrans - based on Fouries series approximation       [MSc Thesis, chapter 9]
    % g. quadgk     - built-in MATLAB Gauss-Kronrod quadrature  [MSc Thesis, chapter 10a]
    % h. hilbert    - built-in MATLAB functions                 [MSc Thesis, chapter 10b]
    
    %% Constructor:
    % wrn - display warnings [default: true]

    %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw, 2011-2012]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform used to 
    % better understand and solve the Kramers-Kronig relations in nonlinear optics"
    % krzysztof.parjaszewski@gmail.com

    properties
        h_transforms;
        show_warning;
        t_gcf;
        t_waitbar;
    end
    
    methods % Constructor
      function obj = HilbertTransform(wrn)
        if nargin < 1, wrn = true; end;
        obj.show_warning = wrn;
        obj.h_transforms = struct( ...
          'name', { ...
            'htran' , ...      % MSc Thesis, chapter  4
            'hncX', ...        % MSc Thesis, chapter  5
            'htrancc', ...     % MSc Thesis, chapter  6
            'hfthilbert', ...  % MSc Thesis, chapter  7
            'herhtrans', ...   % MSc Thesis, chapter  8
            'fourhtrans', ...  % MSc Thesis, chapter  9
            'quadgk', ...      % MSc Thesis, chapter 10a
            'hilbert'}, ...    % MSc Thesis, chapter 10b
          'test_1d', { ...  % One-dimensional tests
            @obj.t_htran_1d, ...
            @obj.t_hncX_1d, ...
            @obj.t_htrancc_1d, ...
            @obj.t_hfthilbert_1d, ...
            @obj.t_herhtrans_1d, ...
            @obj.t_fourhtrans_1d, ...
            @obj.t_quadgk_1d, ...
            @obj.t_hilbert_1d}, ...
          'test_2d', { ... % Two-dimensional tests
            @obj.t_htran_2d, ...
            @obj.t_hncX_2d, ...
            @obj.t_htrancc_2d, ...
            @obj.t_hfthilbert_2d, ...
            @obj.t_herhtrans_2d, ...
            @obj.t_fourhtrans_2d, ...
            @obj.t_quadgk_2d, ...
            @obj.t_hilbert_2d} ...
          );
        if (get(0,'CurrentFigure')==0)
          obj.t_gcf = figure;
        else
          obj.t_gcf = gcf;
        end
        obj.t_waitbar = waitbar(0,'Please wait...');
      end
    end

    methods % One dimensional tests
        function t_htran_1d(obj, model) % Paper by I.J. Weinbert
          % Read model
          [x, X, f, l, model_name] = Utils.prepareTest(model);
          
          % Longer values range
          Y = f(X);               

          % Calculate the Hilbert transform values
          HLi = htran(X, real(Y)); 
          HLr = htran(X, imag(Y));
          
          % Rotate the minus/plus problem for models
          [HLr, HLi] = Utils.rotate(model.orientation, HLr, HLi);
          
          % Shorten previously enlonged range
          [Hi, Hr] = Utils.shortenRange(HLi, HLr, l);
          
          % Show final summary
          obj.make1dSummary(x, f(x), Hi, Hr, 't_htran_1d', model_name);
        end
          
        function t_hncX_1d(obj, model) % Newton-Cotes quadrature 
          % Read model
          [x, ~, f, ~, model_name] = Utils.prepareTest(model); l = length(x);
          
          % hncX can maximally work on model with 50 points:
          lX = min(50, l);
          xX = linspace(min(x), max(x), lX);
          
          % Calculate the Hilbert transform values
          Hi = hncX(@(n) real(f(n)), min(xX), max(xX), 10^(-1), 4, 0.1, lX, obj.show_warning, obj.t_waitbar);
          Hr = hncX(@(n) imag(f(n)), min(xX), max(xX), 10^(-1), 4, 0.1, lX, obj.show_warning, obj.t_waitbar); 
          
          % Rotate the minus/plus problem for models
          [Hr, Hi] = Utils.rotate(model.orientation, Hr, Hi);
          
          % We do not longen-shorten range, because this method takes long

          % Show final summary
          obj.make1dSummary(xX, f(xX), Hi, Hr, 't_hncX_1d', model_name);
        end
           
        function t_htrancc_1d(obj, model) % Hilbert Transform with Clenshaw-Curtis 
          % Read model
          [x, ~, f, ~, model_name] = Utils.prepareTest(model);

          % Calculate the Hilbert transform values
          Hi = htrancc(@(n) real(f(n)), x, 10^(-6), obj.t_waitbar);
          Hr = htrancc(@(n) imag(f(n)), x, 10^(-6), obj.t_waitbar);
          
          % Rotate the minus/plus problem for models
          [Hr, Hi] = Utils.rotate(model.orientation, Hr, Hi);
          
          % We do not longen-shorten range, because this method takes long

          % Show final summary
          obj.make1dSummary(x, f(x), Hi, Hr, 't_htrancc_1d', model_name);
        end
       
        function t_hfthilbert_1d(obj, model) % Hilbert Transform with Hartley transform
           % Read model
          [x, X, f, l, model_name] = Utils.prepareTest(model);
          
          % Longer values range
          x = linspace(min(x),max(x),length(x) * 20);
          X = linspace(min(X),max(X),length(X) * 20);
          l = l * 20;
          Y = f(X);               

          % Calculate the Hilbert transform values 
          HLi = hfthilbert(real(Y)); 
          HLr = hfthilbert(imag(Y)); 
          
          % Rotate the minus/plus problem for models
          [HLr, HLi] = Utils.rotate(model.orientation, HLr, HLi);
          
          % Shorten previously enlonged range
          [Hi, Hr] = Utils.shortenRange(HLi, HLr, l);
          
          % Show final summary
          obj.make1dSummary(x, f(x), Hi, Hr, 't_hfthilbert_1d', model_name);
        end
        
        function t_herhtrans_1d(obj, model) % Discrete Hermite transform and Hermite transform 
          % Read model
          [x, X, f, l, model_name] = Utils.prepareTest(model);
          
          % Shorten arguments range by factor of 100 due to computation error
%           scale = 100;
%           newF = @(t) f(scale*t);
%           newX = linspace(min(X)/scale,max(X)/scale,length(X));
%           newx = linspace(min(x)/scale,max(x)/scale,length(x));

          % Calculate the Hilbert transform values
          HLi = herhtrans(@(n) real(f(n)), X, 10^(-3), 14); 
          HLr = herhtrans(@(n) imag(f(n)), X, 10^(-3), 14); 
          
          % Rotate the minus/plus problem for models
          [HLr, HLi] = Utils.rotate(model.orientation, HLr, HLi);
          
          % Shorten previously enlonged range
          [Hi, Hr] = Utils.shortenRange(HLi, HLr, l);

          % Show final summary
          obj.make1dSummary(x, f(x), Hi, Hr, 't_herhtrans_1d', model_name);
        end
        
        function t_fourhtrans_1d(obj, model) % Fouries series approximation 
           % Read model
          [x, X, f, l, model_name] = Utils.prepareTest(model);

          % Calculate the Hilbert transform values
          HLi = fourhtrans(@(n) real(f(n)), X, 1024, 4096, obj.t_waitbar);
          HLr = fourhtrans(@(n) imag(f(n)), X, 1024, 4096, obj.t_waitbar);
          
          % Rotate the minus/plus problem for models
          [HLr, HLi] = Utils.rotate(model.orientation, HLr, HLi);
          
          % Shorten previously enlonged range
          [Hi, Hr] = Utils.shortenRange(HLi, HLr, l);
          
           % Show final summary
          obj.make1dSummary(x, f(x), Hi, Hr, 't_fourhtrans_1d', model_name);
        end
        
        function t_quadgk_1d(obj, model) % Built-in MATLAB Gauss-Kronrod quadrature 
          % Read model
          [x, X, f, l, model_name] = Utils.prepareTest(model);

          % Calculate the Hilbert transform values
          iter = 0; 
          h = obj.t_waitbar;  
          S = abs(max(x) - min(x))/10000;
          HLi = zeros(1, l); HLr = zeros(1, l);
          waitbar(0, h);
          for v=X
              iter = iter + 1; 
              HLi(iter) =   1 ./ pi .* (quadgk(@(n) real(f(n))./(n-v), -inf, v-S, 'RelTol', 1e-3, 'AbsTol', 1e-3) ...
                                      + quadgk(@(n) real(f(n))./(n-v), v+S, inf, 'RelTol', 1e-3, 'AbsTol', 1e-3));
              HLr(iter) =   1 ./ pi .* (quadgk(@(n) imag(f(n))./(n-v), -inf, v-S, 'RelTol', 1e-3, 'AbsTol', 1e-3) ...
                                      + quadgk(@(n) imag(f(n))./(n-v), v+S, inf, 'RelTol', 1e-3, 'AbsTol', 1e-3));
              waitbar(iter/l, h);
          end; 
          
          % Rotate the minus/plus problem for models
          [HLr, HLi] = Utils.rotate(model.orientation, HLr, HLi);
         
          % Shorten previously enlonged range
          [Hi, Hr] = Utils.shortenRange(HLi, HLr, l);

          % Show final summary
          obj.make1dSummary(x, f(x), Hi, Hr, 't_quadgk_1d', model_name);
        end
        
        function t_hilbert_1d(obj, model) % Built-in MATLAB discrete Hilbert transform
          % Read model
          [x, X, f, l, model_name] = Utils.prepareTest(model);
          
          % Longer values range
          Y = f(X);               

          % Calculate the Hilbert transform values 
          HLi = -imag(hilbert(real(Y))); 
          HLr = -imag(hilbert(imag(Y))); 
          
          % Rotate the minus/plus problem for models
          [HLr, HLi] = Utils.rotate(model.orientation, HLr, HLi);
          
          % Shorten previously enlonged range
          [Hi, Hr] = Utils.shortenRange(HLi, HLr, l);
          
          % Show final summary
          obj.make1dSummary(x, f(x), Hi, Hr, 't_hilbert_1d', model_name);
        end   
    end

    methods % Two dimensional tests
        function t_htran_2d(~, model) % Paper by I.J. Weinbert
        end  
        
        function t_hncX_2d(~, model) % Newton-Cotes quadrature 
        end
          
        function t_htrancc_2d(~, model) % Hilbert Transform with Clenshaw-Curtis 
            % Read Model
            x = model.x; [X1, X2] = meshgrid(x); l = length(x); % Read arguments
            f = model.f; % Read function
            
            % Prepare two waitbars, this process will be long
            inh = waitbar(0, 'Please wait', 'Name', 'Hilbert Transform with Clenshaw-Curtis');
            h1 = waitbar(0, 'Please Wait...', 'Name', 'General 2D htrancc Progress ');sumtime = 0;            

            % 2D - loop
            H = zeros(size(X1));
            for i = 1:1:l, 
                tStart=tic;
                hfun = @(n) f(n, X2(i));H(i, :) = htrancc(hfun, x, 10^(-3), inh)'; 
                tElapsed=toc(tStart);sumtime = sumtime + tElapsed;
                timeleft = sumtime*(l-i)/i; minleft = floor(timeleft/60); secleft = timeleft - minleft*60;
                waitbar(i/l, h1, ['General progress: ', num2str(i/l*100), '%, time left: ', num2str(minleft), 'm ', num2str(secleft), 's' ]);
            end; close(h1);close(inh);
            
            % Present results:
            Hhtrancclog = Utils.doLog(H);          
            figure, mesh(X1, X2, Hhtrancclog);
        end
   
        function t_hfthilbert_2d(~, model) % Hilbert Transform with Hartley transform
        end
        
        function t_herhtrans_2d(~, model) % Discrete Hermite transform and Hermite transform 
        end
        
        function t_fourhtrans_2d(~, model) % Fouries series approximation 
        end
        
        function t_quadgk_2d(~, model) % Built-in MATLAB Gauss-Kronrod quadrature 
        end
        
        function t_hilbert_2d(~, model) % Built-in MATLAB functions 
        end
    end
    
    methods % Utils
      function make1dPlot(obj, x, y, cy, hcy, tit_str)
        clf;
        plot(x, y, x, cy, x, hcy); % plot main value and co-values - analytical and hilbert transformed
        hold on;
        error_plot = cy - hcy; % Absolute error
        plot(x, error_plot, 'linewidth', 4, 'color', 'yellow'); % plot yellow error line
        hold off; 
%         f_name = [datestr(now, 'yyyy-mm-dd HH-MM-SS'), ' ', tit_str]; 
        title(tit_str); % make title
%         saveas(obj.t_gcf,f_name,'jpg');
      end
      
      function make1dSummary(obj,args, vals, Hi, Hr, f_title, model_name)
           obj.make1dPlot(args, real(vals), imag(vals), Hi, [f_title, ' im_part_err ', model_name]);
           obj.make1dPlot(args, imag(vals), real(vals), Hr, [f_title, ' re_part_err ', model_name]);  
      end
    end

end