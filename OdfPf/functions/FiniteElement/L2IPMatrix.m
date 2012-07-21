function l2ip = L2IPMatrix(gqrule, npqp)
% L2IPMatrix - Form matrix for L2 inner product.
%   
%   USAGE:
%
%   l2ip = L2IPMatrix(gqrule, npqp)
%
%   INPUT:
%
%   qrule is a QRuleStructure,
%         it is the global quadrature rule associated with
%         the underlying mesh
%   npqp  is m x n, (sparse) 
%         it is the matrix which takes nodal point values
%         to quadrature point values
%
%   OUTPUT:
%
%   l2ip is n x n, (sparse) 
%        the inner product matrix associated with the mesh
%
ngqp = length(gqrule.wts);
i    = 1:ngqp;
wmat = sparse(i, i, gqrule.wts, ngqp, ngqp);
l2ip = npqp' * wmat * npqp;
%
%  Force symmetry on l2ip, since roundoff may make this asymmetric.
%
l2ip = 0.5 * (l2ip + l2ip');
