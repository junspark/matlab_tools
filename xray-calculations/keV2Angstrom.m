function Angstrom = keV2Angstrom(keV)
% keV2Angstrom - converts x-ray energy in keV to wavelength in angstroms
%
%   USAGE:
%
%   Angstrom = keV2Angstrom(keV)
%
%   INPUT:
%
%   keV
%       x-ray energy in keV
%
%   OUTPUT:
%
%   Angstrom
%       x-ray wavelength in anstroms

hc  = 12.39842; % PER X-RAY BOOKLET 2009
Angstrom    = hc./keV;