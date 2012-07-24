function c = GetSHCoeffs(f, frmesh, V)
% GetSHCoeffs - Find matrix relating ODF to weighted pole figure.
%   
%   USAGE:
%
%   c = GetSHCoeffs(f, frmesh, V)
%
%   INPUT:
%
%   f	  	field over orientation space
%   frmesh	is a MeshStructure on orientation space
%   V		spherical harmonic modes
%
%   OUTPUT:
%
%   c		spherical harmonic coefficients
%
%   see DiscreteHarmonics
c   = f*frmesh.l2ip*V;
