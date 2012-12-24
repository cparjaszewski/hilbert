function makehtrancctest
    model = struct('x',linspace(-20,20,256),'f',@(x) ((1+x.*1i)./(x.^2+1)));
    
    % Read model
    [x, X, f, l] = prepareTest(model);

    % Calculate the Hilbert transform values
    HLi = htrancc(@(n) real(f(n)), X, 10^(-6));
    HLr = htrancc(@(n) imag(f(n)), X, 10^(-6));

    [Hi, Hr] = shortenRange(HLi, HLr, l);

    % Show final summary
    make1dSummary(x, f(x), Hi, Hr, 't_hilbert_1d');

end

 function [Hi, Hr] = shortenRange(HLi,HLr,l)
    Hi = HLi(l*2/5:l*3/5-1);
    Hr = HLr(l*2/5:l*3/5-1);
  end

function [x, X, f, l] = prepareTest(model)
      x = model.x; % Read arguments
      f = model.f; % Read function

      % Prepare longer arguments range
      a = min(x); b = max(x);
      d = b - a;
      A = a - 2 * d; B = b + 2 * d;
      l = 5*length(x); 
      X = linspace(A, B, l); % Longer arguments range
  end

function make1dPlot(x, y, cy, hcy, tit_str)
      figure; % create a new figure
      plot(x, y, x, cy, x, hcy); % plot main value and co-values - analytical and hilbert transformed
      hold on; 
      error_plot = cy - hcy;
      plot(x, error_plot, 'linewidth', 4, 'color', 'yellow'); % plot yellow error line
      hold off; 
      title(tit_str); % make title
end
      
function make1dSummary(args, vals, Hi, Hr, f_title)
   make1dPlot(args, real(vals), imag(vals), -Hi, [f_title, ' - error for imag part']);
   make1dPlot(args, imag(vals), real(vals), Hr, [f_title, ' - error for real part']);  
end