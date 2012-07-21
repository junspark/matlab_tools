function rmesh = RefineMesh(mesh, n, sym)
% RefineMesh - Refine a simplex-based mesh.
%   
%   USAGE:
%
%   rmesh = RefineMesh(mesh, n)
%   rmesh = RefineMesh(mesh, n, sym)
%
%   INPUT:
%
%   mesh is a MeshStructure,
%        using 1D, 2D, or 3D simplicial elements
%   n    is a scalar, (positive integer)
%        the subdivision parameter for each simplex
%   sym  is 4 x s, (optional) 
%        the array of quaternions giving the symmetry
%        group for meshes on orientation space; this argument 
%        is not needed if the mesh has a 'symmetries' component
%
%   OUTPUT:
%
%   rmesh is a MeshStructure,
%         it is derived from the input mesh by subdividing each 
%         simplicial element into a number of subelements 
%         based on the subdivision parameter n
%  
%   NOTES:
%
%   *  This routine only subdivides the simplex.  If the 
%      mesh points are mapped, as in sphere meshes, then
%      that mapping needs to be applied after this routine.
%   
if (nargin < 3)
  sym = [];
  if isfield(mesh, 'symmetries')
    sym = mesh.symmetries;
  end
end
%
crd = mesh.crd;
dc  = size(crd, 1);
con = mesh.con;
[d1, ne] = size(con);
dim = d1 - 1; % dimensionality of simplex
%
%  Subdivide reference simplex.
%
[refcrd, refcon] = SubdivideSimplex(dim, n);
%
%  Generate points.
%
pts = SpreadRefPts(mesh, refcrd); % contains duplicates
pts = reshape(pts, [dc ne*size(refcrd, 2)]);
%
%  Create connectivity.
%
[mrc, nrc] = size(refcon);
maxrc = max(refcon(:));  % = size(refcrd, 2);
nums = (0:(ne-1))*maxrc; % node number offsets
nums = repmat(nums, [mrc*nrc, 1]);
nums = reshape(nums, [mrc nrc*ne]);
crep = repmat(refcon, [1 ne]);
cond = crep + nums;      % connectivity with duplicates
%
%  Remove duplicates.
%
[upts, ord, iord] = UniqueVectors(pts);
con = iord(cond);
%
%  Create MeshStructure.
%
rmesh = MeshStructure(upts, con);
if (~isempty(sym))
  rmesh = ReduceMesh(rmesh, sym, 1.0e-7);
  rmesh.symmetries = sym;
end
rmesh.ord = ord;
%


