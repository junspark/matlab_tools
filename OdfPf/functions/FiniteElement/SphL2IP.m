function l2ip = SphL2IP(mesh, qrule)
% SphL2IP - L^2 inner product matrix for sphere.
%   
%   l2ip = SphL2IP(mesh, qrule)
%
%   mesh is a MeshStructure,
%	a mesh on the sphere in any dimension
%   qrule is a QRuleStructure,
%	a quadrature rule on the reference element
%
%   l2ip is n x n, (sparse)
%	it gives the matrix of the L^2 inner product
%	for functions on the mesh
%
%   Notes:
%
%   *  This is essentially SphGQRule reduced to return only
%      the inner product matrix
%
ne        = size(mesh.con, 2);
[dim, nq] = size(qrule.pts);
npts      = nq*ne;
%
gij  = MetricGij(SphDifferential(mesh, qrule.pts));
jac  = sqrt(DetMatArray(gij));
%
wts0 = reshape(repmat(qrule.wts, [1 ne]), [1 npts]);
wts  = jac .* wts0;
%
npqp = NpQpMatrix(mesh, qrule);
i    = 1:npts;
wmat = sparse(i, i, wts, npts, npts);

% form matrix
l2ip = npqp' * wmat * npqp;

% symmetrize due to roundoff
l2ip = 0.5*(l2ip + l2ip'); 
