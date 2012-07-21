function pb = PolyBasis(dim, deg)
% POLYBASIS - Multi-indices for homogeneous polynomial basis.
%
%   VERSION:  $Id: PolyBasis.m 170 2010-03-01 00:29:25Z boyce $
%   
%   STATUS:  in development
%
%   USAGE:
%
%   pb = PolyBasis(dim, deg)
%
%   INPUT:
%
%   dim is a scalar integer,
%	the dimension of the underlying Euclidean space
%   deg is a scalar integer,
%       the degree of the polynomial space
%
%   OUTPUT:
%
%   pb is dim x n,
%       a list of multi-indices (list of exponents) which span the space of
%	polynomials of exactly degree 'deg'; n is the number of basis
%       elements, depending on dimension and degree
%
%   NOTES:
%
if (dim == 1)
  pb = deg;
  return
end
%
%  Higher dimensions.
%
pb = [];
for d=deg:-1:0
  pbd = PolyBasis(dim - 1, d);
  n   = size(pbd, 2);
  pb  = [pb, [pbd; repmat(deg-d, [1 n])]];
end
%
