function noneqv = VerifyEqv(mesh, sym, tol)
% VerifyEqv - Verify the equivalence array for a mesh.
%   
%   USAGE:
%
%   noneqv = VerifyEqv(mesh, sym)
%   noneqv = VerifyEqv(mesh, sym, tol)
%
%   INPUT:
%
%   mesh is a MeshStructure,
%        on Rodrigues parameters
%   sym  is 4 x n, 
%        a list of quaternions for the symmetry group
%   tol  is a scalar, (optional, default = 5.0e-8) 
%        specifying the tolerance for comparison of equivalent 
%        orientations
%
%   OUTPUT:
%
%   noneqv is 1 x l, 
%          a list of indices in the `mesh.eqv' array
%          which fail the verification
%
%   NOTES:
%
%   *  This computes the misorientation angles between the pairs of
%      nodes specified in the equivalence list.  If any angles 
%      exceed the specified tolerance, the corresponding indices
%      in the mesh.eqv array are returned.
%
crd = mesh.crd;
eqv = mesh.eqv;
%
if (isempty(eqv))
  noneqv = [];
  return
end
%
if (nargin < 3) % use default tolerance
  tol = 5.0e-8;
end
%
q1 = QuatOfRod(crd(:, eqv(1, :)));
q2 = QuatOfRod(crd(:, eqv(2, :)));
%
a = Misorientation(q1, q2, sym);
%
noneqv = find(a > tol);
