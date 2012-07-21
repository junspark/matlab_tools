function h1sip = H1SIPMatrixSph(mesh, qrule, grad, varargin)
% H1SIPMATRIX - Form the matrix operator for H1 semi-norm over a
%   MeshStructure on a unit sphere.
%   
% USAGE
%     1) h1sip = H1SIPMatrixSph(mesh, qrule, grad)
%     2) h1sip = H1SIPMatrixSph(mesh, qrule, grad, gij)
%
% INPUTS
%     1) mesh is a MeshStructure on the unit sphere in dimension d
%        containing l elements, m nodes and n nodal equivalences.
%     2) qrule is a QRuleStructure for a reference (elemental) quadrature
%        rule containing p points.
%     3) grad is an NpQpGradMatrix corresponding to the mesh and qrule --
%        this represents the scalar gradient operator and must be of
%        dimensions (d * l * p) x (m - n).
%     4) gij is the (sparse) matrix which takes nodal point values
%        to quadrature point values
%
% OUTPUTS
%     1) h1sip is (m - n) x (m - n), (sparse) the gradient inner-product
%        matrix associated with the mesh.
%
%   SEE ALSO:  NpQpGradMatrix, MetricGij
%
con  = mesh.con;			% connectivity
crd  = mesh.crd;			% coord's (quaternions)

dim  = size(mesh.crd, 1) - 1;		% manifold dimension
%
gqrule = SphGQRule(mesh, qrule);
gwts  = gqrule.wts;
%
nel = size(con, 2);
nqp = length(qrule.wts);
ngqp = nqp*nel;
%
if nargin < 4
    refpts = qrule.pts;
    tngVec = SphDifferential(mesh, refpts);
    gij = MetricGij(tngVec);
elseif nargin == 4
    gij = varargin{1};
end
%
gwts = reshape(repmat(gwts, [dim*dim, 1]), [dim dim ngqp]);
%
wmat = SparseOfMatArray(gij.*gwts);
%
h1sip = grad'*wmat*grad;
% symmetrize...
h1sip = 0.5*(h1sip + h1sip');
