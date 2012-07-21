function m = CubBaseMesh()
%   CubBaseMesh - Return base mesh for cubic symmetries
%  
%   USAGE:
%
%   m = CubBaseMesh
%
%   INPUT:  no inputs
%
%   OUTPUT:
%
%   m is a MeshStructure,
%        on the cubic fundamental region
%
m = LoadMesh('cub-base');
m.symmetries = CubSymmetries;
