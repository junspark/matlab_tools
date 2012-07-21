function diff = RodDifferential(rmesh, refpts)
% RodDifferential - Differential map for Rodrigues
%   
%   USAGE:
%
%   diff = RodDifferential(rmesh, refpts)
%
%   INPUT:
%
%   mesh   is a mesh,
%          on a Rodrigues mesh
%   refpts is 4 x n,
%          a list of points in the reference element,
%          usually the quadrature points, given in barycentric
%          coordinates
%
%   OUTPUT:
%
%   diff is 4 x 3 x nq,
%        a list of tangent vectors at each reference
%        point for each element; nq is the number of global
%        quadrature points, that is n x ne, where ne is the
%        number of elements
%
nref = size(refpts, 2);
nele = size(rmesh.con, 2);
%
pts  = reshape(SpreadRefPts(rmesh, refpts), [3, nref*nele]);
%
all123 = rmesh.crd(:, rmesh.con(1:3, :));
all4   = rmesh.crd(:, rmesh.con(4, :));
alltan = all123 - reshape(repmat(all4, [3 1]), [3 3*nele]);
alltan = reshape(alltan, [9 nele]);
%
d1 = reshape(repmat(alltan, [nref, 1]), [3 3 nref*nele]);
%
d2 = RDiff(pts);
%
diff = MultMatArray(d2, d1);
%
%-------------------- Local Functions
%
function rd = RDiff(rod)
% RDIFF - Differential of Rodrigues to Quaternion
%
%  Not accounting for map ref-element to mesh.
%   
nrod = size(rod, 2);
rd = zeros(4, 3, nrod);
%
cos_th2 = 1./sqrt(1 + dot(rod, rod));
%
tmp1 = reshape([ones(1, nrod); rod], [4 1 nrod]);
tmp2 = permute(reshape(rod, [3 1 nrod]), [2 1 3]);

term1 = reshape(repmat(-cos_th2.^3, [12 1]), [4 3 nrod]).* ...
	MultMatArray(tmp1, tmp2);
%
term2 = reshape(repmat(cos_th2, [12 1]), [4 3 nrod]).* ...
	repmat([0 0 0; eye(3)],[1 1 nrod]);
%
rd = term1 + term2;
