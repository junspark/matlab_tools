function fun = PolyBasisEval(pbasis, coeffs, pts)
% PolyBasisEval - Evaluate homogeneous polynomials.
%
%   VERSION:  $Id: PolyBasisEval.m 170 2010-03-01 00:29:25Z boyce $
%
%   STATUS:  in development
%
%   USAGE:
%
%   fun = PolyBasisEval(pbasis, coeffs, pts)
%
%   INPUT:
%
%   pbasis is d x m
%          a list of multi-indices representing a
%          polynomial basis for some dimension and degree;
%          d is the dimension of the space; m is the number of 
%          basis elements for the given degree
%
%   coeffs is an m-vector
%          it represents a linear combination of the
%          the basis vectors
%
%   pts    is d x k
%          a list of k points of dimension d
%
%   OUTPUT:
%
%   fun is 1 x k
%       a list of values of the polynomial at the k points
%
%   NOTES:
%
npts = size(pts, 2);
fun = zeros(1, npts);
nba = size(pbasis, 2);
for p=1:npts
  pt  = repmat(pts(:, p), [1 nba]).^pbasis;
  fun(p) = dot(coeffs, prod(pt, 1));
end
