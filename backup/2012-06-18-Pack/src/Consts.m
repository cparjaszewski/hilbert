classdef Consts < handle
    properties
      ia; % accuracy of integration over radius
      ee; % number of exponent elements
      zr; % [mm] range of the z-scan
      d0; % [mm] distance between the z=0 and aperture plane
      ra; % [mm] aperture radius
      T;  % T-factor
      k;  % [1/mm] wave vector
      nop;% number of points
    end
    
    methods
        function obj = Consts(lambda, nop)
           obj.ia = 100;             
           obj.ee = 5;               
           obj.zr = 60;              
           obj.d0 = 600;             
           obj.ra = 1.0;             
           obj.T = 0.0; 
           obj.k = 2*10^6*pi/lambda;
           obj.nop = nop;
        end
    end
end