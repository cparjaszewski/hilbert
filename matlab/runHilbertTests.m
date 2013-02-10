function runHilbertTests
  %% Description:
  % runHilbertTests - tests runner

  %% Test info:
  % In this test in a double all defined models from the AllModels class 
  % are cross-tested with all hilbert transform functions from the the
  % HilbertTransform class.
  
  %% Tested methods:
  % 'htran' 
  % 'hncX 
  % 'htrancc' 
  % 'frthilbert' 
  % 'herhtrans'
  % 'fourhtrans' 
  % 'quadgk' 
  % 'hilbert' 

  %% Author info:
  % [Krzysztof Parjaszewski, University of Wroclaw]
  % As a part of MSc Thesis - "Numerical evaluation of the Hilbert transform used to 
  % better understand and solve the Kramers-Kronig relations in nonlinear optics"
  % krzysztof.parjaszewski@gmail.com
  
  % Initialization
  clear;
  
  % Starting up
  ht = HilbertTransform(false); 
  allmodels = AllModels();
  % We iterate through all hilbert transform' methods:
  for htransform = ht.h_transforms, 
    tStart=tic; % timing
    % We iterate through all models:
    for model = allmodels.all_models, 
     if strcmp('htrancc', htransform.name), htransform.test_1d(model); end
%      htransform.test_1d(model); 
    end; 
    tElapsed=toc(tStart); % timing
    testSummaryText = ['Test for ', htransform.name, ' took: ', num2str(tElapsed), '.'] ;
    display(testSummaryText);
  end;
  
  % Ending up
  close(ht.t_waitbar);
  close(gcf);
end