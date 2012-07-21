function mesh = MeshStructure(crd, con, eqv)
% MESHSTRUCTURE - Create mesh structure from mesh data.
%   
%   mesh = MeshStructure(crd, con, eqv)
%   mesh = MeshStructure(crd, con)
%   mesh = MeshStructure
%
%   crd is e x n, the nodal point array 
%   con is d x m, the connectivity
%   eqv is 2 x k, (optional) the equivalence array
%
%   mesh is a mesh-structure, which consists of the fields
%        above; with no args, it returns and empty mesh
%        structure; with two args, it sets the equivalence 
%        array to be empty
%
if (nargin == 0)
  crd = [];
  con = [];
  eqv = [];
elseif (nargin == 2)
  eqv = [];
end
%
mesh = struct('crd', crd, 'con', con, 'eqv', eqv);
%
