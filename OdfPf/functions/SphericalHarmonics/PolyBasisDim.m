function dspace = PolyBasisDim(dim, deg)
% PolyBasisDim - return dimension of polynomial space
%
%   VERSION:  $Id: PolyBasisDim.m 170 2010-03-01 00:29:25Z boyce $
%
%   STATUS:  in development
%
%   USAGE:
%
%   dspace = PolyBasisDim(dim, deg)
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
%   dspace is an integer
%          the dimension of the space of polynomials of degree 'deg'
%          in dimension 'dim'
%
%   NOTES:
%
%   
if deg == 0
  dspace = 1;
  return
elseif deg < 0
  dspace = 0;
  return
end


switch dim
 case 1
  dspace = 1;
 case 2
  dspace = deg + 1;
 case 3
  dspace = ((deg+1)*(deg+2))/2;
 case 4
  dspace = ((deg+1)*(deg+2)*(deg+3))/6;
 case 5
  dspace = ((deg+1)*(deg+2)*(deg+3)*(deg+4))/24;
 otherwise
  error('dimension out of range')
end
