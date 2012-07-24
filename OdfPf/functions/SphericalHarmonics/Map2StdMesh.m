function [fStd] = Map2StdMesh(f, fMesh, StdMesh)
% Map2StdMesh - projects a field in orientation space from one mesh to another mesh
%   
%   USAGE:
%
%   fStd = Map2StdMesh(f, fMesh, StdMesh)
%
%   INPUT:
%
%   f	  	field over original orientation space mesh
%   fMesh	is the original MeshStructure of orientation space
%   StdMesh is the target MeshStructure of orientation space
%
%   OUTPUT:
%
%   fStd	field over target orientation space mesh
[fele, fcrd]    = FiberCoordinates(StdMesh.crd(:,1:StdMesh.numind), fMesh);
fStd            = EvalMeshFunc(fMesh, f, fele, fcrd);