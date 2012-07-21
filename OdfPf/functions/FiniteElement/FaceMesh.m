function fmesh = FaceMesh(mesh)
% FACEMESH - Compute boundary mesh.
%   
%   fmesh = FaceMesh(mesh)
%
%   mesh is a MeshStructure
%
%   fmesh is a MeshStructure, having the same
%            nodal points and equivalences, but
%            with a boundary connectivity 
%
%   FaceMesh works by extracting element faces and
%   using only those with multiplicity one.
%
fmesh = mesh;
[f, m] = MeshFaces(mesh.con);
fmesh.con = f(:, m == 1);
