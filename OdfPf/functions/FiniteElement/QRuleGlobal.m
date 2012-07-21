function gqrule = QRuleGlobal(mesh, qrule, MetricFun)
% QRuleGlobal - Construct a global quadrature rule for a mesh.
%   
%   USAGE:
%
%   gqrule = QRuleGlobal(mesh, qrule)
%   gqrule = QRuleGlobal(mesh, qrule, @MetricFun)
%
%   INPUT:
%
%   mesh  is a MeshStructure
%   qrule is a QRuleStructure, 
%         the quadrature rule for the reference element 
%  
%   MetricFun is a function handle, (optional)
%             It returns the volumetric integration factors
%             for the applied mapping.  It is for the case
%             in which the mesh is the parameter space for
%             some manifold and has a subsequent mapping
%             associated with each point.
%            
%             met = Metricfun(pts)
%
%             pts is m x n array of points
%             met is 1 x n vector of integration factors
%
%   OUTPUT:
%
%   NOTES:
%
%   *  For spheres, this routine is superceded by `SphGQRule',
%      which handles the metric dependency on the mesh.  This
%      still is to be used for meshes on Rodrigues parameters,
%      for which the metric function is independent of the mesh,
%      or for meshes with no further mapping involved.
%             
crd = mesh.crd; ndim = size(crd, 1);
con = mesh.con; nele = size(con, 2);
%
qpt = qrule.pts; nqpt = size(qpt, 2);
qwt = qrule.wts;
%
ngqp = nqpt*nele;
gqp  = reshape(SpreadRefPts(mesh, qpt), [ndim ngqp]);
%
%  Weights consist of Jacobian * local weights * Metric.
%
lwts = repmat(qwt, [1 nele]);
ejac = reshape(repmat(Jacobian(mesh), [nqpt  1]), [1 ngqp]);
%
if (nargin > 2)
  gwt = lwts .* ejac .* feval(MetricFun, gqp);
else
  gwt = lwts .* ejac;
end
%
gqrule = QRuleStructure(gqp, gwt);
