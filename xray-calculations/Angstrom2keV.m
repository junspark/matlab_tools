function keV = Angstrom2keV(Angstrom)
% Angstrom2keV - converts x-ray wavelength in angstrom to energy in keV
%
%   USAGE:
%
%   keV = Angstrom2keV(Angstrom)
%
%   INPUT:
%
%   Angstrom
%       x-ray wavelength in angstroms
%
%   OUTPUT:
%
%   keV
%       x-ray energy in keV

hc  = 12.39842; % PER X-RAY BOOKLET 2009
keV = hc./Angstrom;