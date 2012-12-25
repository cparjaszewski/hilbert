function [x,y,z] = plot3dCoumarine307(wavelengths, immagGamma)
% getting the wavelengths:
minwl = wavelength(1)/1.4; 
maxwl= wavelength(2)*1.4;

% prepare arguments for integration
args = minwl:0.1E-08:maxwl;

% prepare the mesh grid
[XX,YY] = meshgrid(args,args);






