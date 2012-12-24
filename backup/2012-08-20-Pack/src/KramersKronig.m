classdef KramersKronig  < handle 
% KramersKronig v.1.0.0
% Matlab class to enable solving linear and nonlinear kramers kronig
% relations
%
% Application written as a part of "Organometallics in Nanophotonics" project
% and my MSc thesis
% 
% K.Parjaszewski 
% http://www.organometallics.pwr.wroc.pl/krzysztof.parjaszewski
% University of Wroclaw, Institute of Computer Science
% 15th November 2010

    properties % Variables
        complexFun;
    end
     
    methods % Constructor
        function obj = KramersKronig(cFun)
           obj.complexFun = cFun;
        end 
    end
    
    methods % Calculations
        function res = runCPV(obj)
            res = obj.complexFun(1);
        end
        
        function res = runQuadkGK(obj)
            res = obj.complexFun(1);
        end
        
        function res = runCC(obj)
            res = obj.complexFun(1);
        end
        
        function res = runRK(obj)
            res = obj.complexFun(1);
        end
        
        function res = runKeller(obj)
           res = obj.complexFun(1);
        end
    end
    
end