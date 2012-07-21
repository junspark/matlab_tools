function [gqrule, l2ip] = SphGQRule(mesh, qrule)
% SphGQRule - Global quadrature rule for sphere.
%   
%   USAGE:
%
%   [gqrule, l2ip] = SphGQRule(mesh, qrule)
%
%   INPUT:
%
%   mesh  is a MeshStructure,
%         on a sphere of any dimension
%   qrule is a QRuleStructure,
%         on the reference element of the mesh
%
%   OUTPUT:
%
%   gqrule is a QRuleStructure,
%          for the entire mesh
%   l2ip   is n x n, (sparse)
%          it gives the l2 inner product in terms of function
%          values at the nodal points
%
ne = size(mesh.con, 2);
[dim, nq] = size(qrule.pts);
npts = nq*ne;
%
diff = SphDifferential(mesh, qrule.pts);
gij  = MetricGij(diff);
jac  = sqrt(DetMatArray(gij));
%
wts0 = reshape(repmat(qrule.wts, [1 ne]), [1 npts]);
%
pts = reshape(SpreadRefPts(mesh, qrule.pts), [dim npts]);
pts = UnitVector(pts);
wts = jac .* wts0;
%
gqrule = QRuleStructure(pts, wts);
%
if (nargout >= 2)
  npqp = NpQpMatrix(mesh, qrule);
  i = 1:npts;
  wmat = sparse(i, i, wts, npts, npts);
  l2ip = npqp' * wmat * npqp;
end
