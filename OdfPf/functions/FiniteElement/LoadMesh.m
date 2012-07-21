function mesh = LoadMesh(name)
% LoadMesh - Load mesh from a ascii file.
%  
%   USAGE:
%
%   mesh = LoadMesh(name)
%
%   INPUT:
%
%   name is a string, 
%       the basename of the mesh files
%
%   OUTPUT:
%
%   mesh is a MeshStructure, 
%        for the mesh being loaded
%
%   NOTES:
%
%   * Expected file suffixes are:
%
%        .crd for nodal point coordinates,
%        .con for connectivity
%        .eqv for equivalences (optional)
%   
try
  crd = load ([name '.crd'])';
  con = load ([name '.con'])';
catch
  error('failed to open mesh data files');
end
%
try
  eqv = load ([name '.eqv'])';
catch
  eqv = [];
end
%
mesh = MeshStructure(crd, con, eqv);
