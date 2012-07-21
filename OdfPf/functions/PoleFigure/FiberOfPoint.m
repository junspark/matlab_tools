function fib = FiberOfPoint(p, h, ndiv, qsym, invfib)
% FiberOfPoint - Find fiber above point for specified pole figure.
%   
%   USAGE:
%
%   fib = FiberOfPoint(p, h, ndiv, qsym)
%   fib = FiberOfPoint(p, h, ndiv, qsym, invfib)
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
%   invfib is a scalar, (optional) 
%        a flag which, if nonzero, produces the inverse fiber 
%   
%   OUTPUT:
%
%   fib is 3 x ndiv x n, 
%       the column represents the Rodrigues
%       vector of a point on the fiber, the row spans the points
%       on the fiber, and the page spans the points on the sphere
%   
%   NOTES:
%
%   * `h' need not be a unit vector; it is normalized here
%
if (nargin < 5) 
  invfib = 0;
end
%
n = size(p, 2);
ntot = n*ndiv;
%
qfib = QFiberOfPoint(p, h, ndiv, qsym, invfib);
qfib = reshape(qfib, 4, ntot);
%
fib = ToFundamentalRegion(qfib, qsym);
fib = reshape(fib, 3, ndiv, n);
