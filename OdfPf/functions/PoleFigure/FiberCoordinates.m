function [fele, fcrd] = FiberCoordinates(fib, mesh)
% FiberCoordinates - Find mesh coordinates for a given fiber.
%   
%   USAGE:
%
%   [fele, fcrd] = FiberCoordinates(fib, mesh)
%
%   INPUT:
%
%   fib is a 3 x ndiv x npole,
%       it is an array with each page giving the fiber over 
%       the n'th pole point; usually the result of FiberOfPoint
%   mesh is a MeshStructure 
%       on the fundamental region
%
%   OUTPUT:
%
%   fele is ndiv x npole, 
%        the list of elements containing the fiber points.
%   fcrd is 4 x ndiv x npole, 
%        the barycentric coordinates of the fiber points
%
%  Notes:
%
%  *  The call to the matlab builtin `tsearchn' can fail, possibly
%     due to the fact that the meshes are not necessarily Delaunay
%     tesselations.
%
crd = mesh.crd;
con = mesh.con;
%
npole = size(fib, 3);
ndiv  = size(fib, 2);
ndim  = size(fib, 1);% should be three
npts  = npole*ndiv;
%
%  This section requires some error checking since 
%  tsearchn may fail to produce the barycentric coordinates.
%  
fibpts = reshape(fib, ndim, npts)';
[tmpele, tmpcrd] = tsearchn(crd', con', fibpts);
%
failed = ~isfinite(tmpele);
if (any(failed))       % Retry.  This can work because a different
                       % algorithm is used for large point sets.
  [tmpele(failed), tmpcrd(failed, :)] = ...
      tsearchn(crd', con', fibpts(failed, :));
  if (any(~isfinite(tmpele)))
    error('Failed to find fiber coordinates.')
  end
end
%
fele = reshape(tmpele, ndiv, npole);
fcrd = reshape(tmpcrd', [4 ndiv npole]);
