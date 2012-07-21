function p2d = EqualArea(p3d, basis)
% EQUALAREA - Equal area projection of 3D points on sphere.
%   
%   p2d = EqualArea(p3d)
%   p2d = EqualArea(p3d, basis)
%
%   p3d   is a real 3 x n array:
%            a list of n unit vectors
%   basis is a real 3 x 3 matrix:  (optional)
%            it's columns form an orthonormal basis used to find the
%            projection; defaults to the identity
%
%   p2d is a real 2 x n array:
%          the equal area projections of points of `p3d'
%
%   The equal area projection is computed using the 
%   third basis vector as the pole.  Planar components
%   are given relative to the first two basis vectors.
%
if (nargin == 1)
  pcrd = p3d;
else
  pcrd = basis'*p3d;     % components in basis
end
%
nrmxy = dot(pcrd(1:2, :), pcrd(1:2, :));      
nonzr = (nrmxy > eps);
rscal = nrmxy + (1 - pcrd(3, :)).^2;          
rscal (nonzr)  = rscal(nonzr)./nrmxy(nonzr);   
rscal (~nonzr) = 1.0;
rscal = sqrt(rscal);
%
p2d = repmat(rscal, [2 1]) .* pcrd(1:2, :);
