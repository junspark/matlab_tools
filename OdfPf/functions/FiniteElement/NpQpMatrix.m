function npqp = NpQpMatrix(mesh, qrule)
% NpQpMatrix - Nodal point values to quadrature points.
%   
%   USAGE:
%
%   npqp = NpQpMatrix(mesh, qrule)
%
%   INPUT:
%
%   mesh  is a MeshStructure
%   qrule is a QRuleStructure, 
%         a quadrature rule on the reference element
%
%   OUTPUT:
%
%   npqp is m x n, (sparse) 
%        it is the matrix which takes the values
%        at the independent nodes of the mesh to the values 
%        at the quadrature points;  here, n is the number of 
%        independent nodal values, and m = q*e, where e is the 
%        number of elements in the mesh and q is the number of 
%        qudrature points per element
%                  
crd = mesh.crd;
con = mesh.con;
eqv = mesh.eqv;
%
qpt = qrule.pts;
qwt = qrule.wts;
%
ncrd = size(crd, 2); 
scon = size(con)   ; nnpe = scon(1); nele = scon(2); 
nqpt = size(qpt, 2); % size(qpt, 1) should match nnpe
%
ngqp = nele*nqpt;
%
%npqp = sparse(ngqp, nred)
nqij = repmat(reshape(qpt, [1 nnpe*nqpt]), [1 nele]);
nq_i = reshape(repmat((1:ngqp), [nnpe 1]), [1 ngqp*nnpe]);
nq_j = reshape(repmat(con, [nqpt 1]), [1 ngqp*nnpe]);
%
npqp = sparse(nq_i, nq_j, nqij, ngqp, ncrd); % needs to be reduced
%
if (~isempty(eqv))
  npqp = EqvReduce(npqp, eqv); % reduce along rows
end
