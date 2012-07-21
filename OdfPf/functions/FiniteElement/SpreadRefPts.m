function pts = SpreadRefPts(mesh, ref)
% SpreadRefPts - Spread reference coordinates to all elements.
%   
%   USAGE:
%
%   pts = SpreadRefPts(mesh, ref)
%
%   INPUT:
%
%   mesh is a MeshStructure 
%        with simplicial element type
%   ref is e x k, 
%       a list of barycentric coordinates for points in
%       the reference element; 
%
%   OUTPUT:
%
%   pts is d x k x m, 
%       the spatial coordinates of the reference points under
%       the isoparametric mapping defined by the mesh;
%       `d' is the dimension of the mapped space, 'k' is the
%       number of reference points, and 'm' is the number
%       of elements in `mesh'
%
%   NOTES:
%
%   *  Typically this would be used for numerical quadrature.
%   *  e = d + 1 is not necessary; cross-dimension mappings are ok.
%
crd = mesh.crd;
con = mesh.con;
%
d = size(crd, 1); % dimension of points
%e = size(con, 1); % nodes per element 
%
nref = size(ref, 2); % number of reference points
if (nref == 0)
  pts = [];
  return
end
%
nels = size(con, 2); % number of elements
%
pts = zeros(d, nref, nels);
for i=1:nels;
  pts(:, :, i) = crd(:, con(:,i)')*ref;
end
%
