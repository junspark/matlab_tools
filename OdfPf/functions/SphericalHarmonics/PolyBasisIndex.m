function ind = PolyBasisIndex(mi)
% PolyBasisIndex - Find index of multi-index in given basis.
%
%   VERSION:  $Id: PolyBasisIndex.m 170 2010-03-01 00:29:25Z boyce $
%
%   STATUS:  in development
%
%   USAGE:
%
%   ind = PolyBasisIndex(mi)
%
%   INPUT:
%
%   mi is n x k,
%	a list of k multi-indices in n-dimensions; they
%	are expected to be all of the same degree
%
%   OUTPUT:
%
%   ind is 1 x k,
%       the list of basis indices
%
%   NOTES:
%
%   *  Currently, this only works up to dimension n=4.
%
[n, k] = size(mi);  % dimension of space
if k == 0
  ind = [];
  return
end
%
deg = sum(mi, 1);
%
for j=1:k
  indj = 0;
  mydeg = cumsum(mi(:, j));
  for d=n:-1:2
    indj = indj + PolyBasisDim(d, mydeg(d)) ...
	        - PolyBasisDim(d, mydeg(d) - mi(d, j) );
  end
  ind(j) = indj + 1;
end

