% Inputs to the code
% ______________________________
% Input the radial locations of diffraction measurements (mm) in an
% - in increasing order
R   = 6.35 + [
    0.15; 
    0.29; 
    0.43; 
    0.58; 
    0.74; 
    0.98; 
    1.2; 
    1.5; 
    1.9; 
    2.5; 
    3.3; 
    4.4];


% Input the coordinates along the thickness (mm)
% - in increasing order
t   = [-0.58 -0.48 -0.38 -0.28 -0.18 -0.08 0.02 0.12 0.22 0.32 0.42 0.52];
dt  = 0.05;
%
% Alfa measurements (deg)
% - in increasing order
% - Elements will be symmetrically positioned around that angle
% - if there is single angle please enter its variation
alpha=[90];
% Variation if there is a single alpha (not necessary for multiple alphas)
dalpha=5;