function qfib = QFiberOfPoint(p, h, ndiv, qsym, invfib)
% QFiberOfPoint - Find fiber above point in quaternions.
%   
%   USAGE:
%
%   qfib = QFiberOfPoint(p, h, ndiv, qsym)
%   qfib = QFiberOfPoint(p, h, ndiv, qsym, invfib)
%   
%   INPUT:
%
%   p    is 3 x n,  
%        an array of points on the sphere 
%   h    is 3 x 1,  
%        a specified pole direction; if `invfib' is nonzero,
%        this is interpreted as a crystal direction; otherwise
%        it is interpreted as a sample direction
%   ndiv is a positive integer, 
%        the number of equally spaced divisions along the fiber
%   qsym is 4 x m, 
%        the quaternion representation of the symmetry group
%   invfib is a scalar, (optional, default: 0) 
%        a flag which, if nonzero, produces the inverse fiber 
%   
%   OUTPUT:
%
%   qfib is 4 x ndiv x n, 
%        each column gives the quaternion representation of a point 
%        on the fiber; the row spans the points on the fiber;
%        and the page (third index) spans the points on the sphere
%   
%   NOTES:
%      
%   * h need not be a unit vector; it is normalized here
%
if (nargin < 5)
  invfib = 0;
end
%
n = size(p, 2); % number of points
cutoff = 1.0e-8;
%
%  First, normalize h to make it a unit vector.
%
h = h./norm(h);
% 
%  To find the fiber, find one orientation above p, then
%  premultiply by any rotation about h.
%
hn = repmat(h, [1 n]);
%
%  We start with a rotation of pi about the midpoint.
%
ax = p + hn;
anrm = sqrt(dot(ax, ax, 1));
okay = (anrm > cutoff);
nokay = sum(okay);
if (nokay == n)
  ax = ax./repmat(anrm, 3, 1);
else
  %
  %  Need to find a vector orthogonal to h.   
  %
  %hcross = cross(repmat(h, [1 3]), eye(3));
  %[nmax, imax] = max(dot(hcross, hcross, 1));
  %hperp = hcross(:, imax)/sqrt(nmax);
  nspace = null(h');
  hperp  = nspace(:, 1);
  %
  %   Make sure the arrays are not empty before
  %   handling the two cases separately.
  %
  if nokay > 0
    ax(:, okay)  = ax(:, okay)./repmat(anrm(okay), 3, 1);
  end
  %
  %  We know there are some, not OK.
  %
  ax(:, ~okay) = repmat(hperp, 1, n-nokay); 
  %
end
%
q0 = [zeros(1, n); ax];
%
%  Now find rotations which leave h fixed. 
%
phi = (0:(ndiv-1))*((2*pi)/ndiv);
qh  = QuatOfAngleAxis(phi, repmat(h, [1 ndiv]));
%
qfib = zeros(4, ndiv, n);
for i=1:n
  qfib(:, :, i) = QuatProd(repmat(q0(:, i), [1 ndiv]), qh);
end
%
if (invfib ~= 0) % inverse fiber
  qfib(2:4, :, :) = -qfib(2:4, :, :);
end
