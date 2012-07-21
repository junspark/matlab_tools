function plp = PolyLaPlacian(mi)
% POLYLAPLACIAN - LaPlacian of polynomial (multi-indices)
%   
%   VERSION:  $Id: PolyLaPlacian.m 170 2010-03-01 00:29:25Z boyce $
%
%   STATUS:  in development
%
%   USAGE:
%
%   plp = PolyLaPlacian(mi)
%
%   INPUT:
%
%   mi is d x n
%      an array of multi-indices of dimension d (should 
%      be a full space of multi-indices)
%
%   OUTPUT:
%
%   plp is d x n (sparse)
%       an array of coefficients for the Polynomial Basis of 
%       degree two less than that of the input (mi);
%
%   NOTES:
%
%   * It may be more useful to return a full array.
%   * Needs fixing for the case when all inputs have zero LaPlacian
%
[n, k]   = size(mi);
deg   = sum(mi(:, 1));
lpdim = PolyBasisDim(n, deg - 2);
%
mi_ind = 1:k;
plp = sparse(lpdim, k);
%
for d=1:n
  mi_d = mi(d, :);
  coef = mi_d .* (mi_d - 1);
  lmid = mi;
  lmid(d, :) = mi(d,:) - 2;
  %
  nonz = (coef ~= 0);
  indd = PolyBasisIndex(lmid(:, nonz));
  plp  = plp + sparse(indd, mi_ind(nonz), coef(nonz), lpdim, k);
end
