classdef Utils < handle
    %% Description:
    % Utils - a model with common functions widely used in 

    %% Author info:
    % [Krzysztof Parjaszewski, University of Wroclaw]
    % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform in~nonlinear optics"
    % krzysztof.parjaszewski@gmail.com
  
   methods (Static)
      % Removes nans and infinites
      function res = getres(condition, val)
          res = condition .* val;
          res = Utils.onlynums(res);
      end
      
      % Zeros NaN and infinities
      function res = onlynums(input)
          res = input;
          res(isnan(res) | (isfinite(res)==0)) = 0;
      end
      
      % We are rotating the imaginary/real part of result
      function [Rout, Iout] = rotate(orientation, Rin, Iin)
        
        switch(orientation)
          case 1,
            Rout = - Rin;
            Iout = Iin;
          case 2,
            Rout = Rin;
            Iout = - Iin;
          case 3,
            Rout = Rin;
            Iout = Iin;
          case 4,
            Rout = Rin;
            Iout = Iin;
          otherwise
            Rout = Rin;
            Iout = Iin;
        end
        
      end
      
      % Make the logarithmic scale for positive and negative values
      function output = doLog(input)
          output  = zeros(size(input));
          upper   = input; upper(input < 0) = 0;
          lower   = input; lower(input > 0) = 0; lower = - lower;
          output  = output + log(upper+1);
          output  = output - log(lower+1);
          output  = Utils.getres(output);        
      end
      
      % Common test initialization for one-dimensional model
      function [x, X, f, l, model_name] = prepareTest(model)
          x = model.x; % Read arguments
          f = model.f; % Read function
          model_name = model.name;

          % Prepare longer arguments range
          a = min(x); b = max(x);
          d = b - a;
          A = a - 2 * d; B = b + 2 * d;
          l = 5*length(x); 
          X = linspace(A, B, l); % Longer arguments range
      end
      
      % Shorten previously enlonged range
      function [Hi, Hr] = shortenRange(HLi,HLr,l)
        Hi = HLi(l*2/5:l*3/5-1);
        Hr = HLr(l*2/5:l*3/5-1);
      end
          
    end
end