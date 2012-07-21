function m = HexBaseMesh()
%   HexBaseMesh - Return base mesh for hexagonal symmetries
%  
%   USAGE:
%
%   m = HexBaseMesh
%
%   INPUT:  none
%
%   OUTPUT:
%
%   m is a MeshStructure,
%     on the hexagonal fundamental region
%
m = LoadMesh('hex-base');
m.symmetries = HexSymmetries;
