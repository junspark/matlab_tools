function m = OrtBaseMesh()
%   OrtBaseMesh - Base mesh for orthorhombic symmetry.
%  
%   USAGE:
%
%   m = OrtBaseMesh
%
%   INPUT:  no inputs
%
%   OUTPUT:
%
%   m is a MeshStructure,
%        on the orthorhombic fundamental region
%
m = LoadMesh('ort-base');
