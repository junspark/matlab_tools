function CheckMesh(mesh, vector)
% CheckMesh - Check that array sizes match for mesh structure.
%   
%   USAGE:
%
%   CheckMesh(mesh)
%   CheckMesh(mesh, vector)
%
%   INPUT:
%
%   mesh   is a MeshStructure
%   vector is a vector, (optional) 
%          a reduced vector of values on the mesh nodal points
%
%   OUTPUT: calls `error' if inconsistency is found
%
%   NOTES:
%
%   This routine checks that the coordinates and equivalence
%   array of a mesh are consistent.  If the vector argument
%   is given, it checks that it is of the correct length.
%   If any array sizes are inconsistent, it prints an error
%   message and returns.
%
crd = mesh.crd;
eqv = mesh.eqv;
%
ncrd = size(crd, 2);
if (isempty(eqv))
  maxeqv = ncrd;
else
  maxeqv = max(eqv(1, :));
end
%
if (ncrd ~= maxeqv)
  error('Mesh coordinates inconsistent with equivalence array.')
end
%
if (nargin == 2)
  if (isempty(eqv))
    nreduced = ncrd;
  else
    nreduced = min(eqv(1, :)) - 1;
  end
  lenvec = length(vector);
  if (lenvec ~= nreduced)
    error ('Vector does not match mesh.');
  end
end
