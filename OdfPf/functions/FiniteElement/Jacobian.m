function jac = Jacobian(mesh)
% Jacobian - Compute Jacobian of linear mesh mappings.
%
%   USAGE:
%
%   jac = Jacobian(mesh)
%
%   INPUT:
%
%   mesh is a MeshStructure,
%        with simplicial element type
%
%   OUTPUT:
%
%   jac is 1 x m, 
%       the Jacobian of each element
%
%   NOTES:
%
%   *  The mesh may be embedded in a space of higher 
%      dimension than the reference element.  In that
%      case, the Jacobian is computed as (sqrt(det(J'*J))
%      and is always positive.  When the target space is
%      of the same dimension as the reference element,
%      the Jacobian is computed as usual and can be
%      positive or negative.
%
%   *  Only simplicial (linear) element types are allowed.
%
crd = mesh.crd;
con = mesh.con;
%
e    = size(con, 1); % number of nodes/element
ddom = e - 1;        % dimension of domain space
dtar = size(crd, 1); % number of dimensions in target space
%
nels = size(con, 2); % number of elements
%
jac = zeros(1, nels);
%
if (ddom == dtar)  % map to same dimensional space
  for i=1:nels
    simp = crd(:, con(:, i)');
    mat  = simp(:, 1:ddom) - repmat(simp(:, e), 1, ddom);
    jac(i) = det(mat);
  end
else               % embedding in higher dimensions
  for i=1:nels
    simp = crd(:, con(:, i)');
    mat  = simp(:, 1:ddom) - repmat(simp(:, e), 1, ddom);
    mat  = mat'*mat;
    jac(i) = sqrt(det(mat));
  end
end

