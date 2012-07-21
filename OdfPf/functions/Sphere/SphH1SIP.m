function h1sip = SphH1SIP(mesh, qrule)
% SphH1SIP -- H^1 semi-inner product on sphere
%
%   USAGE:
%
%   h1sip = SphH1SIP(mesh, qrule)
%
%   INPUT:
%
%   mesh  is a MeshStructure,
%         a mesh on a sphere of any dimension
%   qrule is a QRuleStructure,
%         a quadrature rule for the reference simplex
%
%   OUTPUT:
%
%   h1sip is n x n, (sparse)
%         it is the matrix of the H^1 semi-inner product,
%         where n is the number of independent degrees of
%         freedom associated with the mesh
%        
refpts = qrule.pts;
gqrule = SphGQRule(mesh, qrule);
%
dim    = size(mesh.con, 1) - 1;
dim1   = dim + 1;
nrps   = size(refpts, 2);
nels   = size(mesh.con, 2);
ncrd   = size(mesh.crd, 2);
neqv   = size(mesh.eqv, 2);
%
refder = RefDerivatives(dim);
gij    = MetricGij(SphDifferential(mesh, refpts));
npts   = size(gij, 3);
refder = repmat(refder, [1 1 npts]);
dercmp = MultMatArray(InvMatArray(gij), refder);
%
%  Set up indices.
%
ndnum = repmat(mesh.con, [nrps 1]);
ndnum = ndnum(:)';
ndnum = repmat(ndnum, [dim 1]);
ndnum = ndnum(:);
%
qpnum = 1:nrps*nels;
qpnum = repmat(qpnum, [dim1 1]);
qpnum = qpnum(:)';
qpdof1 = dim * (qpnum - 1) + 1;
qpdof = qpdof1;
for d = 2:dim
  qpdof = [qpdof; qpdof1 + d - 1];
end
qpdof = qpdof(:);
%
gnpqp = sparse(qpdof, ndnum, dercmp(:), dim*npts, ncrd);
if (~isempty(mesh.eqv))
  gnpqp = EqvReduce(gnpqp, mesh.eqv);
end
%  Now form semi-inner product matrix.
%
wts = reshape(repmat(gqrule.wts, [dim*dim, 1]), [dim dim npts]);
gij = SparseOfMatArray(gij.*wts);
h1sip = gnpqp' * gij * gnpqp;
%
%  The inner product matrix should be symmetric,
%  but it can lose symmetry due to roundoff.  It 
%  needs to be symmetric for eigenvalue computations.
%
h1sip = 0.5 * (h1sip + h1sip');
%%
%if (nargout > 1)
%  grad = gnpqp;
%end
