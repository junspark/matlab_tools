% Inputs to the code
% ______________________________
% Input the radial locations of diffraction measurements (mm) in an
% - in increasing order
R   = 6.35 + ...
    [.15 .285 .427 .58 .745 .92 1.11 1.315 1.54 1.78 2.06 2.36 2.7 3.07 3.5 4 4.6 5.3 6.2 7.25 8.7]';

% Input the coordinates along the thickness (mm)
% - in increasing order
% t   = [-0.5 -0.25 0 0.25 0.5];
t   = linspace(-0.5, 0.5, 2);
dt  = 0.05;
%
% Alfa measurements (deg)
% - in increasing order
% - Elements will be symmetrically positioned around that angle
% - if there is single angle please enter its variation
alpha=[180 157.5 135 112.5 90];

% Variation if there is a single alpha (not necessary for multiple alphas)
dalpha=5;
return% Inputs to the code
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
t   = [-0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5];
dt  = 0.05;
%
% Alfa measurements (deg)
% - in increasing order
% - Elements will be symmetrically positioned around that angle
% - if there is single angle please enter its variation
alpha=[180 157.5 135 112.5 90];
alpha=[90];
% Variation if there is a single alpha (not necessary for multiple alphas)
dalpha=5;
return
