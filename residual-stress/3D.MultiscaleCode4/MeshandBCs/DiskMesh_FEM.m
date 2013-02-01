% This code is written to generate the mesh and the boundary conditions for
% the shrink-fit-disk problem
clear
close all
clc
%
% Inputs (diffraction data centers)
inputs
%
% Calculate DV stresses using the analytical solution
[ndv,dvc_x,dvc_y,dvc_z]=findnodecoordinates(R,t,alpha);
%
%
%
% Values of R, t, and alpha at the nodes
[Rnodes,tnodes,alphanodes] = findnodalvalues(R,t,dt,alpha,dalpha);
%
%
% Find the coordinates of the nodes
[numnp,x,y,z]=findnodecoordinates(Rnodes,tnodes,alphanodes);
%
%
% Find the connectivity matrix
[np] = findconnectivity(Rnodes,tnodes,alphanodes);
numel=size(np,1);
% 
% Apply boundary conditions
[nps,numels,npsyms,numelsyms]=applyBCs_FEM(Rnodes,tnodes,alphanodes,x,y,z);
%
%
% Plot the mesh
%plotmesh
%
%
% Mesh properties (8-node brick)
nnps=4;
meltyp=8;
numeq = 3*numnp;
nnpe = meltyp;
%
%
% Save the meshdata
save meshdataFEM numnp numeq numel numels meltyp nnpe nnps np nps x y z
%    dvc_x dvc_y dvc_z ndv
%
%