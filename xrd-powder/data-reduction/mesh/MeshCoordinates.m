function [els, crds] = MeshCoordinates(mesh, pts)
% MeshCoordinates - find elements and barycentric coordinates of points
%   
%   USAGE:
%
%   [els, crds] = MeshCoordinates(mesh, pts)
%
%   INPUT:
%
%   mesh is a MeshStructure,
%        it should be Delaunay, but it may be okay anyway if it is
%        not too irregular
%   pts  is m x n,
%        a list of n m-dimensional points
%
%   OUTPUT:
%
%   els is an n-vector, (integers)
%       the list of element numbers containing each point
%   crds is (m+1) x n,  
%       the barycentric coordinates of each point
%
%   NOTES:
%
%   *  Only available for simplicial element types.
%
%   *  The call to the matlab builtin `tsearchn' can fail if 
%      the mesh is not Delaunay.
%
MYNAME = 'MeshCoordinates';
%
crd = mesh.crd;
con = mesh.con;
%
[m, n] = size(pts);
%
%  This section requires some error checking since 
%  tsearchn may fail to produce the barycentric coordinates.
%  
[tmpele, tmpcrd] = tsearchn(crd', con', pts');
%
failed = isnan(tmpele);
%
if (any(failed))       % Retry.  This can work because a different
                       % algorithm is used for large point sets.
  warning([MYNAME, ':SearchFailed'], ...
	  'Search failed for %d points.  Retrying.', sum(failed))
  [tmpele(failed), tmpcrd(failed, :)] = ...
      tsearchn(crd', con', pts(:, failed)');
  if (any(isnan(tmpele)))
    error('Failed to find fiber coordinates.')
  end
end
%
els  = tmpele;
%
%  This is to handle matlab inconsistency in tsearchn.
%
if (n == 1) 
  crds = tmpcrd;
else
  crds = tmpcrd';
end
