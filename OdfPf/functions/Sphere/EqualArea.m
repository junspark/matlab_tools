function p2d = EqualArea(p3d, basis)
% EqualArea - Equal area projection on sphere.
%   
%   USAGE:
%
%   p2d = EqualArea(p3d)
%   p2d = EqualArea(p3d, basis)
%
%   INPUT:
%
%   p3d   is 3 x n,
%         a list of n unit vectors
%   basis is a 3 x 3 matrix, (optional)
%         it's columns form an orthonormal basis used to find the
%         projection; defaults to the identity
%
%   OUTPUT:
%
%   p2d is a real 2 x n array:
%       the equal area projections of points in `p3d'
%
%   NOTES:
%
%   *  The equal area projection is computed using the 
%      third basis vector as the pole.  Planar components
%      are given relative to the first two basis vectors.
%
if (nargin == 1)
  pcrd = p3d;
else
  pcrd = basis'*p3d;     % components in basis
end
%
zfac = sqrt(2./(1 + pcrd(3, :)));
p2d  = repmat(zfac, [2, 1]) .* pcrd(1:2, :);
