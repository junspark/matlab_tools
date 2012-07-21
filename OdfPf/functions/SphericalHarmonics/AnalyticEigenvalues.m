function [eVal, mult] = AnalyticEigenvalues(dim, deg)
% AnalyticEigenvalues - Analytic spherical eigenvalues.
%   
%   VERSION:  $Id: AnalyticEigenvalues.m 170 2010-03-01 00:29:25Z boyce $
%
%   STATUS:  in development
%
%   USAGE:
%
%   eVal = AnalyticEigenvalues(dim, deg)
%
%   INPUT:
%
%   dim is an integer
%          the dimension of the sphere 
%
%   deg is an integer
%          the degree of the polynomial space
%
%   OUTPUT:
%
%   eVal is a scalar
%           the eigenvalue of the space
%   mult is an integer
%           the multiplicity of eigenvectors (dimension of eigenspace)
%
%   NOTES:
%
%
switch dim
 case 1
  eVal = deg^2;
  mult = 2;
 case 2
  eVal = deg*deg + deg;
  mult = 2*deg + 1;
 case 3
  eVal = deg^2 + 2*deg;
  mult = (deg + 1)^2;  
 otherwise
  error('dimension out of range:  needs to be 3 or less')
end
