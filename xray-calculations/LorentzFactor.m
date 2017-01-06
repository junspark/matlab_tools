function L = LorentzFactor(th)
% LorentzFactor - Lorentz factor calculation
%
%   USAGE:
%
%   L = LorentzFactor(th)
%
%   INPUT:
%
%   th
%       angle in deg
%
%   OUTPUT:
%
%   L
%       lorentz factor
sth     = sind(th);
s2th    = sind(2.*th);

L   = 1./sth;
L   = L./s2th;
