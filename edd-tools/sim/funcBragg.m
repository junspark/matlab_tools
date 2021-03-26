function lambda = funcBragg(tth, d_hkl)
% funcBragg - Bragg function to determine the 2 theta angle in EDD setup
%
%   USAGE:
%
%   y = funcBragg(tth, d_hkl)
%
%   INPUT:
%
%   tth
%       take off angle (degrees)
% 
%   d_hkl
%       known d-spacings (Angstroms)
%
%   OUTPUT:
%
%   lambda
%       wavelength computed based on tth and d_hkl (Angstrom)
%
lambda  = 2 * d_hkl * sind (tth/2);