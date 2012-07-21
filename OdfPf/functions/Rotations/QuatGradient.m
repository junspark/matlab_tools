function g = QuatGradient(m, qrule)
% QuatGradient - Build gradient on quaternion manifold
%   
%   VERSION:  $Id$
%
%   STATUS:  in development
%
%   USAGE:
%
%   g = QuatGradient(m, qrule)
%
%   INPUT:
%
%   m     is a MeshStructure
%            a mesh on Rodrigues space
%
%   qrule is a QuadratureRule
%
%
%   OUTPUT:
%
%   g is 4 x nsf x nqp x nel
%        the gradient of each shape function at each QP of each element
%
%   NOTES:
%
dim   = 3; 
numEl = size(m.con, 2);
numQp = length(qrule.wts);
%
%  Find derivative components
%
tVecs    = RodDifferential(m, qrule.pts);
gij      = MetricGij(tVecs);
refDeriv = RefDerivatives(dim);
%
nEqp = size(gij, 3);
g = zeros(4, 4, nEqp);
for eqp=1:nEqp
  g(:, :, eqp) = tVecs(:, :, eqp)*(gij(:, :, eqp)\refDeriv);
end
%
g = reshape(g, [4 4 numQp numEl]);
