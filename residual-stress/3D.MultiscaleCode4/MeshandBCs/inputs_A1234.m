% Inputs to the code
% ______________________________
% Input the radial locations of diffraction measurements (mm) in an
% - in increasing order
R   = 6.35 + [
    0.500;
    0.900;
    1.450;
    1.950;
    2.400;
    2.950];


% Input the coordinates along the thickness (mm)
% - in increasing order
t   = [-0.525 -0.375 -0.225 -0.075 0.075 0.225 0.375 0.525];
dt  = 0.05;
%
% Alfa measurements (deg)
% - in increasing order
% - Elements will be symmetrically positioned around that angle
% - if there is single angle please enter its variation
alpha=[90 120 150 180];
% Variation if there is a single alpha (not necessary for multiple alphas)
dalpha=5;
