function h1sip = H1SIPMatrix(mesh, qrule, diff)
% H1SIPMatrix -- H^1 semi-inner product [Documentation out of date!]
%
%   USAGE:
%
%   h1sip = H1SIPMatrix(mesh, qrule, diff)
%
%   INPUT:
%
%   mesh  is a MeshStructure
%         the mesh
%   qrule is a QRuleStructure
%         a qudrature rule (local)
%   diff is d1 x d2 x n
%        the list of tangent vectors at each of n points
%
%   OUTPUT:
%
%   h1sip is n x n, sparse
%       it is the matrix of the H^1 semi-inner product,
%	where n is the number of independent degrees of
%	freedom associated with the mesh
%
%   NOTES:
%
%   * Modified on 2009-01-20 to improve memory usage; this routine
%     was very limiting for large meshes; 
%   --- clearing variables when no longer needed
%        
[dimT, dim, npts] = size(diff);
dim1 = dim + 1;
nele = size(mesh.con, 2);
ncrd = size(mesh.crd, 2);
neqv = size(mesh.eqv, 2);
nrps = length(qrule.wts);
%
refder = RefDerivatives(dim);
gij    = MetricGij(diff);
refder = repmat(refder, [1 1 npts]);
gqrwts = sqrt(DetMatArray(gij)).* repmat(qrule.wts(:)', [1 nele]);
%
%  Set up indices.
%
ndnum = repmat(mesh.con, [nrps 1]);
ndnum = ndnum(:)';
ndnum = repmat(ndnum, [dim 1]);
ndnum = ndnum(:);
%
qpnum = 1:nrps*nele;
qpnum = repmat(qpnum, [dim1 1]);
qpnum = qpnum(:)';
qpdof1 = dim * (qpnum - 1) + 1;
clear qpnum
qpdof = qpdof1;
for d = 2:dim
  qpdof = [qpdof; qpdof1 + d - 1];
end
clear qpdof1 
qpdof = qpdof(:);
%
%  Generate NP -> QP matrix
%
dercmp = MultMatArray(InvMatArray(gij), refder);
clear refder
gnpqp = sparse(qpdof, ndnum, dercmp(:), dim*npts, ncrd);
clear dercmp
%
if (~isempty(mesh.eqv))
  gnpqp = EqvReduce(gnpqp, mesh.eqv);
end
%  Now form semi-inner product matrix.
%
%wts = reshape(repmat(gqrwts, [dim*dim, 1]), [dim dim npts]);
%gij = SparseOfMatArray(gij.*wts);
gij = SparseOfMatArray(gij.* reshape(repmat(gqrwts, [dim*dim, 1]), [dim dim npts]) );
h1sip = gnpqp' * gij * gnpqp;
%
%  The inner product matrix should be symmetric,
%  but it can lose symmetry due to roundoff.  It 
%  needs to be symmetric for eigenvalue computations.
%
h1sip = 0.5 * (h1sip + h1sip');
%
