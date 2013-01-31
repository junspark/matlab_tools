function [faces, mult] = MeshFaces(con)
% MeshFaces - Find lower dimensional faces for a mesh.
%   
%   *** This has been superseded by MeshBoundary.  ***
%
%   USAGE:
%
%   [faces, (opt) mult] = MeshFaces(con)
%
%   INPUT:
%
%   con is d x n, 
%       the mesh connectivity (tetrahedra or triangles)
%
%   OUTPUT:
%
%   faces is (d-1) x m, (integer)
%         the connectivity of the lower dimensional
%                       faces of con
%   mult is 1 x m,  (integer, optional)
%        the multiplicity of each face; this should
%        be either 1 (for a surface face) or 2 (for
%        an interior face)
%
d = size(con, 1); % dimension plus one
n = size(con, 2); % number of elements
dface = d - 1;    % dimension of face
nallf = d*n;      % total number of faces, including interior ones
%
%  efaces is the connectivity for faces within an element
%
efaces = CycleIndices(d);
efaces = efaces(1:dface, :);
%
f = [];
for i=1:d;
  f = [f    con(efaces(:, i), :)];
end
%
allfaces = sortrows(sort(f)')';
[u, f_of, a_of] = unique(allfaces', 'rows');
%
%  Evaluate output args.
%
faces = u';
%
%  Multiplicity if requested.  
%
%  Note:  can get rid of for-loop as in SymEquiv
%         (or use sparse matrix to accumulate)
%
if (nargout >= 2)
  nfaces = size(faces, 2);
  wunz = ones(1, nallf);
  %
  mult = sparse(wunz, a_of, wunz, 1, nfaces);
  mult = full(mult);   % note sparse accumulates entries
  %
%  for i=1:nallf
%    mult(a_of(i)) = mult(a_of(i)) + 1;
%  end
end
