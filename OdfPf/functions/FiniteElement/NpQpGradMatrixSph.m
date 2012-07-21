function [grad, varargout] = NpQpGradMatrixSph(mesh, refpts)
% NPQPGRADMATRIX - Form the gradient operator for scalar functions 
%                   over the given mesh.
%   
%          grad = NpQpGradMatrixSph(mesh, refpts)
%   [grad, gij] = NpQpGradMatrixSph(mesh, refpts)
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
con  = mesh.con;			% connectivity
crd  = mesh.crd;			% coord's (quaternions)
eqv  = mesh.eqv;			% equivalnce array

dim  = size(mesh.crd, 1) - 1;		% manifold dimension
dim1 = dim + 1;			            % embedded space dim

inv = 1;                    % inverse metric flag

nel  = size(con, 2);			% # elements
nqp  = size(refpts, 2);		    % # quadrature pts/element
ngqp = nqp*nel;			        % # quad pts/mesh

% Form and tile the shape function derivatives
dN_dxi = SimplexSFunDer(dim);
dN_dxi = repmat(dN_dxi, [1, 1, nqp*nel]);

% Form tangent vectors under mapping at each quad pt
tngVec = SphDifferential(mesh, refpts);

% Form the metric tensor and its inverse at every quad pt
[gij, invGij] = MetricGij(tngVec, inv);

% Form array of gradient components (along tangent vectors)
gradCompMatT = MultMatArray(dN_dxi, invGij);
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
