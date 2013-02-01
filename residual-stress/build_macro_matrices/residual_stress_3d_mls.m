% The input mesh and connectivity will be provided by an input file 
% Input .txt file has the format: (is read by mesh_data.m)
% a. total number of elements in the mesh
% b. total number of nodes in the mesh
% c. type of element used in the mesh
% d. coordinates of the nodes
% e. connectivity matrix (ordering of the nodes is important,
% - see element_library.pdf)
% f. nodes of the free surfaces (NEW)
%
clear all
close all
clc
%
% Input constants for the analysis

% 10 is firm constraint
% Penalty parameter for equilibrium
lambda  = 1;
% Penalty parameter for traction-free surfaces
beta    = 0.01;
% Penalty parameter for symmetry surfaces
alpha   = 0.01;
% % Radius for the spline interpolation used in the moving least square (mls)
% % method to interpolate from the stresses at the quad points from the given
% % diffraction measuremeents 
%rho=8.5;
% Scaling factor for the ellipsoids
kappa   = 6;

% Number of diffraction volumes
% diffraction volume grids
% THIS IS MAT FILE FROM STEP3
load /home/jspark/Documents/work/Data/APS201103/GE/A1234.DV.GRID.mat

% Load mesh data (mesh+BCs)
% THIS IS MAT FILE FROM STEP4
load meshdata.mat
% load meshdata_KM.mat
plot3(x,y,z, 'b.')
hold on
plot3(dvc_x, dvc_y, dvc_z, 'k.')

% Calculate the ellipsoids
[evals,evecs]=findevs(dvc_x,dvc_y,dvc_z,ndv);

% evaluates elemental information: shape function derivatives etc at the
% quadrature points
[nqptv,wtq,sfac,dndxi,dndet,dndze,nqpts,swt,ssfac,dnda,dndb] = shafac(meltyp,nnpe,nnps);
%
% Calculate quad point coordinates
[xqpt,yqpt,zqpt] = qptloc(x,y,z,np,numel,nnpe,sfac);
%
% Find the nodal stresses which match the given (or measured) stresses
% while satisfying equilibrium and zero surface traction conditions!
Ksr=zeros((numel*36*nnpe*nnpe)+(36*numelsyms*nnps*nnps)+(36*numels*nnps*nnps),1);
Ksc=zeros((numel*36*nnpe*nnpe)+(36*numelsyms*nnps*nnps)+(36*numels*nnps*nnps),1);
Ksv=zeros((numel*36*nnpe*nnpe)+(36*numelsyms*nnps*nnps)+(36*numels*nnps*nnps),1);

% For mls method
Qsr=zeros(numel*nnpe*36*ndv,1); 
Qsc=zeros(numel*nnpe*36*ndv,1); 
Qsv=zeros(numel*nnpe*36*ndv,1);
fc=1; kc=1; qc=1;

% loop over the elements to set up matrices
for   iele =1:1:numel
    [dndx,dndy,dndz,detj] = sfder(iele,nnpe,nqptv,sfac,dndxi,dndet,dndze,np,x,y,z);

    % Find the weighthing value for the quad points using mls method
    [phi] = mlsel3d(xqpt,yqpt,zqpt,dvc_x,dvc_y,dvc_z,evals,evecs,kappa,iele,nqptv,ndv);
    
    % Residual calculation using mls method 
    [se qe] = elstif_residual_mls(nnpe,nqptv,ndv,wtq,sfac,dndx,dndy,dndz,detj,lambda,phi);

    % Residual calculation without using mls method
    % [se qe] = elstif_residual(nnpe,nqptv,wtq,sfac,dndx,dndy,dndz,detj,lambda);
    [Ksr,Ksc,Ksv,kc] = assmbl_ful(iele,nnpe,np,se,Ksr,Ksc,Ksv,kc);
    
    % Assemble Q matrix (mls method)
    [Qsr,Qsc,Qsv,qc] = assmbl_ful_q_mls(iele,nnpe,ndv,np,qe,Qsr,Qsc,Qsv,qc);
    
    % Assemble Q matrix (without mls method)
    % [Qsr,Qsc,Qsv,qc] = assmbl_ful_q(iele,nnpe,np,qe,Qsr,Qsc,Qsv,qc);

    disp(['L2 fit + eq. constraint: ' num2str(iele/numel*100) '% is complete!'])
end
%

% T MATRIX
% free surface condition
% Calculate the constant matrices
[bigNsurf]  = bigNsurfmat(nnps,nqpts,ssfac);
% For each prescribed free surface
for   ieles = 1:1:numels
%
    % Surface normal and jacobian
    [n rjs] = surfjac(ieles,nnps,nqpts,dnda,dndb,nps,x,y,z);
%
    % Form the surface "stiffness matrix"
    sek = freesurf_residual(nnps,nqpts,swt,bigNsurf,n,rjs,beta);
%
    % Assemble the free surface condition to the corresponding location
    % in the overall stiffness matrix
    [Ksr,Ksc,Ksv,kc] = assmblsurf_ful(ieles,nnps,nps,sek,Ksr,Ksc,Ksv,kc);
%
    disp(['Applying free surface constraint: ' num2str(ieles/numels*100) '% is complete!']);
end

% D MATRIX
% symmetry surface condition (in-plane surface tractions are free)
% For each prescribed free surface
for   ieles = 1:1:numelsyms
%
    % Surface normal and jacobian
    [n rjs] = surfjac(ieles,nnps,nqpts,dnda,dndb,npsyms,x,y,z);
%    
    % Form the surface "stiffness matrix"
    sek = symsurf_residual(nnps,nqpts,swt,bigNsurf,n,rjs,alpha);
%
    % Assemble the free surface condition to the corresponding location
    % in the overall stiffness matrix
    [Ksr,Ksc,Ksv,kc] = assmblsurf_ful(ieles,nnps,npsyms,sek,Ksr,Ksc,Ksv,kc);
%
    disp(['Applying symmetry surface constraint: ' num2str(ieles/numelsyms*100) '% is complete!'])
end
%
%
%
%
% put in sparse global format
%
% Fs=sparse(Fsr,1,Fsv);
% 
K=sparse(Ksr,Ksc,Ksv);
%
Q=sparse(Qsr,Qsc,Qsv);
%
%
%  impose boundary conditions  --  not yet done
%
% [Ks,Fs] = bc_ful(numnp,npbc,bfx,bfy,Ks,Fs);
%
% solve 
%
% gsig=full(K\F);
%
%
% Plot the results
% output
%
%