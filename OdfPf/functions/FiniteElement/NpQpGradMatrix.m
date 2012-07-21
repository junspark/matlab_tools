function [grad, varargout] = NpQpGradMatrix(mesh, refpts)
% NPQPGRADMATRIX - Form the gradient operator for scalar functions 
%                   over the given mesh.
%   
%          grad = NpQpGradMatrix(mesh, refpts)
%   [grad, gij] = NpQpGradMatrix(mesh, refpts)
%
%   mesh:   a MeshStructure on the unit sphere in dimension d containing
%             l elements, m nodes and n nodal equivalences
%   refpts: a d x p array  of barycentric coordinates  (typically 
%             quadrature points) in the reference element at which to 
%             compute the gradient 
%
%   grad:   is ((d - 1) * l * p) x (m - n), (sparse) the gradient
%             operator matrix associated with the mesh
%   gij:    is (d - 1) x (d - 1) x p array of metric tensors at each of
%             the p refpts          
%
%   SEE ALSO:  NpQpMatrix
%
con  = mesh.con;			    % connectivity
crd  = mesh.crd;			    % coord's (quaternions)
eqv  = mesh.eqv;			    % equivalnce array

dim  = size(crd, 1);      % manifold dimension
dim1 = dim + 1;			      % embedded space dim

inv = 1;                  % inverse metric flag

nel  = size(con, 2);		  % # elements
nqp  = size(refpts, 2);		% # quadrature pts/element
ngqp = nqp*nel;			      % # quad pts/mesh

% Form and tile the shape function derivatives (wrt to the parent element
% coord's)
dN_dxi = SimplexSFunDer(dim);
dN_dxi = repmat(dN_dxi, [1, 1, nel]);

% % % Form tangent vectors under mapping at each quad pt
% % tngVec = RodDifferential(mesh, refpts);

% Form the metric tensor and its inverse at every quad pt
refpts = SpreadRefPts(mesh, refpts);
gij    = RodMetricGij(refpts(:, :));
invGij = InvMatArray(gij);
% % [gij, invGij] = MetricGij(tngVec, inv);

% Find the elemental Jacobians
elcrds = reshape(crd(:, con), [dim, dim1, nel]);
Jac = MultMatArray(permute(dN_dxi, [2, 1, 3]), permute(elcrds, [2, 1, 3]));

% Cartesian gradient on the mesh space...
dN_dx = MultMatArray(InvMatArray(Jac), permute(dN_dxi, [2, 1, 3]));
dN_dx = reshape(permute(repmat(dN_dx, [1, 1, 1, nqp]), [2, 1, 4, 3]), [dim1, dim, ngqp]);

% Form array of gradient components (along tangent vectors)
gradCompMatT = MultMatArray(dN_dx, invGij);
gradComps = gradCompMatT(:);

% Form row and column index vectors to build sparse array
rowIndex = 1:dim*ngqp;
rowIndex = repmat(rowIndex, [dim1, 1]);
rowIndex = rowIndex(:);
colIndex = repmat(con, [dim*nqp, 1]);
colIndex = colIndex(:);

% Form the gradient operator
grad = sparse(rowIndex, colIndex, gradComps);

% Reduce along rows using equivalences
if (~isempty(eqv))
  grad = EqvReduce(grad, eqv); 
end

% Output metric array if requested
if nargout == 2
  varargout{1} = gij;
end
