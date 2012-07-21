function [uvec, ord, iord] = UniqueVectors(vec, tol)
% UniqueVectors - Remove near duplicates from a list of vectors.
%   
%   USAGE:
%
%   [uvec, ord, iord] = UniqueVectors(vec)
%   [uvec, ord, iord] = UniqueVectors(vec, tol)
%
%   INPUT:
%
%   vec is d x n, 
%       an array of n d-vectors
%   tol is a scalar, (optional) 
%       the tolerance for comparison; it defaults to 1.0e-14
%
%   OUTPUT:
%
%   uvec is d x m, 
%        the set of unique vectors; two adjacent vectors are considered
%        equal if each component is within the given tolerance
%   ord  is an m-vector, (integer)
%        which relates the input vector to the output vector, 
%        i.e. uvec = vec(:, ord)
%   iord is an n-vector, (integer)
%        which relates the reduced vector to the original vector, 
%        i.e. vec = uvec(:, iord)
%
%   NOTES:
%
%   *  After sorting, only adjacent entries are tested for equality
%      within the tolerance.  For example, if x1 and x2 are within
%      the tolerance, and x2 and x3 are within the tolerance, then 
%      all 3 will be considered the same point, even though x1 and
%      x3 may not be within the tolerance.  Consequently, if you
%      make the tolerance too large, all the points will be
%      considered the same.  Nevertheless, this routine should be 
%      adequate for the its intended application (meshing), where
%      the points fall into well-separated clusters.
%
if (nargin < 2)
 tol = 1.0e-14; 
end
%
[d n] = size(vec);
%
%  Sort each row and assign an integer to each
%  independent value.
%
ivec = zeros(d, n);
%
for row=1:d
  %
  [tmpsrt, tmpord] = sort(vec(row, :));
  %
  %  Compare entries and assign rank.
  %
  tmpcmp = abs(tmpsrt(2:n) - tmpsrt(1:(n-1)));
  indep  = logical([1 (tmpcmp > tol)]);  % independent values
  rowint = cumsum(indep);
  ivec(row, tmpord) = rowint;
  %
end
% 
%  Sort the integer vector uniquely and save ordering.
%
[utmp, ord, iord] = unique(ivec', 'rows');
%
%  Return sorted vector and row vectors of indices.
%
uvec  = vec(:, ord);
ord   = ord';
iord  = iord';
