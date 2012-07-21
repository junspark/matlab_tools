function opm = OdfPfMatrix(hkl, mesh, sym, pts, div, invpf)
% OdfPfMatrix - Find matrix relating ODF to pole figure.
%   
%   USAGE:
%
%   opm = OdfPfMatrix(hkl, mesh, sym, pts, div)
%   opm = OdfPfMatrix(hkl, mesh, sym, pts, div, invpf)
%
%   INPUT:
%
%   hkl   is a 3-vector, 
%         the crystal direction specifying the pole figure
%   mesh  is a MeshStructure,
%         on orientation space
%   sym   is 4 x s, 
%         the symmetry group in quaternions
%   pts   is 3 x p, 
%         a list of p points on the sphere (S^2)
%   div   is a positive integer, 
%         the number of divisions per fiber to use in
%         calculating the fiber integral
%   invpf is a scalar, (optional, default = 0) 
%         a nonzero value flags causes computation of inverse pole 
%         figure matrix in which the `hkl' is interpreted
%         as a fixed sample direction and `pts' is an array
%         of crystal directions
%
%   OUTPUT:
%
%   opm is p x n, (sparse) 
%       the matrix which takes nodal point values on the
%       fundamental region to pole figure or inverse
%       values at the specified points
%
%   NOTES:
%
%   *  opm acts on the "reduced" set of nodes, the set of nodes
%      in which equivalent nodes are combined into a single
%      degree of freedom
%   *  this routine can be very memory intensive; you may need to
%      build the matrix in pieces; if you have a lot of points
%      and a large value of `div', then you should break the points
%      into smaller groups 
%   *  the ODF-PF matrix preserves mean value, so that 
%      an ODF in MUD (multiples of uniform) maps to a pole
%      figure also in MUD.
%
if (nargin < 6)
  invpf = 0;
end
%
crd = mesh.crd;
con = mesh.con;
eqv = mesh.eqv;
%
na = size(crd, 2);      % total number of nodes
n  = min(eqv(1,:)) - 1; % the number of nonequivalent nodes
ne = na - n;            % number of equivalences
p  = size(pts, 2);      % number of pole figure points
%
opm   = sparse(p, na);
%
fib = FiberOfPoint(pts, hkl(:), div, sym, invpf);
%
[fele, fcrd] = FiberCoordinates(fib, mesh);
%
nterms    = 4*div;
ntermstot = nterms*p;
i   = reshape(repmat(1:p, nterms, 1), ntermstot, 1);
j   = reshape(con(:, fele),           ntermstot, 1);
vij = reshape(fcrd,                   ntermstot, 1);
opm = sparse(i, j, vij, p, na);
opm = (1/div)*EqvReduce(opm, eqv);
