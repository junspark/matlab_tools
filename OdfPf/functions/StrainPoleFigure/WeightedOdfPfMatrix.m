function Aopm = WeightedOdfPfMatrix(hkl, mesh, sym, pts, div, odf)
% WeightedOdfPfMatrix - Find matrix relating ODF to weighted pole figure.
%   
%   USAGE:
%
%   Aopm = WeightedOdfPfMatrix(hkl, mesh, sym, pts, div, odf)
%
%   INPUT:
%
%   hkl   is a 3-vector, 
%         the crystal direction specifying the pole figure
%   mesh  is a MeshStructure on orientation space
%   sym   is 4 x s, 
%         the symmetry group in quaternions
%   pts   is 3 x p, 
%         a list of p points on the sphere (S^2)
%   div   is a positive integer, 
%         the number of divisions per fiber to use in
%         calculating the fiber integral
%   odf   is 1 x numind,
%         is the odf over mesh 
%
%   OUTPUT:
%
%   Aopm  is p x n, (sparse) 
%         the matrix which takes nodal point values on the
%         fundamental region to odf weighted values pole figure 
%         at the specified points
%
%   NOTES:
%
%   *  Aopm acts on the "reduced" set of nodes, the set of nodes
%      in which equivalent nodes are combined into a single
%      degree of freedom
%   *  this routine can be very memory intensive; you may need to
%      build the matrix in pieces; if you have a lot of points
%      and a large value of `div', then you should break the points
%      into smaller groups 

%
invpf   = 0;
%
crd = mesh.crd;
con = mesh.con;
eqv = mesh.eqv;
%
na  = size(crd, 2);      % total number of nodes
p   = size(pts, 2);      % number of pole figure points
%
fib = FiberOfPoint(pts, hkl(:), div, sym, invpf);
%
[fele, fcrd] = FiberCoordinates(fib, mesh);
%
Avals = EvalMeshFunc(mesh, odf, fele, fcrd);
%
nterms    = 4*div;
ntermstot = nterms*p;
i   = reshape(repmat(1:p, nterms, 1), ntermstot, 1);
j   = reshape(con(:, fele),           ntermstot, 1);
%
Avij    = reshape(repmat(Avals,4,1).*fcrd, ntermstot, 1);
Aopm    = sparse(i, j, Avij, p, na);
Aopm    = (1/sum(Avals))*EqvReduce(Aopm, eqv);
