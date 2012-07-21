function fun = DiscreteDelta(mesh, l2ip, elem, ecrd, wts)
% DiscreteDelta - Evaluate sum of discrete delta functions.
%
%   USAGE:
%
%   fun = DiscreteDelta(mesh, l2ip, elem, ecrd, wts)
%
%   INPUT:
%
%   mesh is a MeshStructure
%   l2ip is n x n,
%        the inner product matrix for the mesh, where
%        n is the number of independent degrees of freedom
%   elem is an n-vector of integers,
%        the list of elements containing the delta-function points
%   ecrd is m x n,
%        the list of barycentric coordinates of the delta-function points
%   wts  is an n-vector,(optional)
%        it gives the weights of each point; if not present,
%        the weights are set to one
%
%   OUTPUT:
%
%   fun is an n-vector,
%       the nodal point values of the sum of the discrete
%       delta functions associated with the points and weights
%
%   NOTES:
%
%   *  The resulting function is not normalized in any way.
%
ncrd = size(mesh.crd, 2);
%
npts   = length(elem);
[m, n] = size(ecrd);
%
if (n ~= npts)
    error('Array size mismatch:  elem and ecrd')
end
%
if (nargin < 5)
    wts = ones(1, n);
end
%
mn = m*n;
i = mesh.con(:, elem); i = i(:);
j = ones(mn, 1);
w = repmat(wts(:)', [m 1]); w = w(:);
v = ecrd(:).*w;
%
rhs = sparse(i, j, v, ncrd, 1);
rhs = full(rhs);
%
%  Now reduce according to equivalence.
%
rhs = EqvReduce(rhs', mesh.eqv);
%
fun = l2ip\rhs';
