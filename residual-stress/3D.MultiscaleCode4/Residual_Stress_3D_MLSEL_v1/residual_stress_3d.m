clear all
close all
clc
%
% Input constants for the analysis
% Penalty parameter for equilibrium
lambda  = 0.1;
% Penalty parameter for symmetry surfaces
alpha   = 10;
% Penalty parameter for traction-free surfaces
beta    = 1;

lambda  = 0.1;
alpha   = 0;
beta    = 0.1;

% % Radius for the spline interpolation used in the moving least square (mls)
% % method to interpolate from the stresses at the quad points from the given
% % diffraction measuremeents 
% rho=6;
% Number of diffraction volumes
% ndv=210;
%
% % Load diffraction data
% % diffvolcentroids: Diffraction volume coordinates
% load DVcent_full.mat
%
% load DVdata.mat
%
% load DVdata105.mat
%
% read in the information for a global mesh
% load meshdata_KM.mat
load meshdata.mat
%
% evaluates elemental information: shape function derivatives etc at the
% quadrature points
[nqptv,wtq,sfac,dndxi,dndet,dndze,nqpts,swt,ssfac,dnda,dndb] = shafac(meltyp,nnpe,nnps);
%
% Calculate quad point coordinates
[xqpt,yqpt,zqpt] = qptloc(x,y,z,np,numel,nnpe,sfac);
%
% Find the nodal stresses which match the given (or measured) stresses
% while satisfying equilibrium and zero surface traction conditions!
Ksr = zeros((numel*36*nnpe*nnpe)+(36*numelsyms*nnps*nnps)+(36*numels*nnps*nnps),1);
Ksc = zeros((numel*36*nnpe*nnpe)+(36*numelsyms*nnps*nnps)+(36*numels*nnps*nnps),1);
Ksv = zeros((numel*36*nnpe*nnpe)+(36*numelsyms*nnps*nnps)+(36*numels*nnps*nnps),1);

% Without mls
Qsr = zeros(numel*36*nnpe,1);
Qsc = zeros(numel*36*nnpe,1);
Qsv = zeros(numel*36*nnpe,1);

qc      = 1;
kc      = 1;
% loop over the elements to set up matrices
for iele = 1:1:numel
    [dndx,dndy,dndz,detj]   = sfder(iele,nnpe,nqptv,sfac,dndxi,dndet,dndze,np,x,y,z);
    
    % Resiudal calculation without using mls method
    [se qe] = elstif_residual(nnpe,nqptv,wtq,sfac,dndx,dndy,dndz,detj,lambda);
    
    % Assemble K matrix (without mls method)
    [Ksr,Ksc,Ksv,kc]                = assmbl_ful(iele,nnpe,np,se,Ksr,Ksc,Ksv,kc);
    
    % Assemble Q matrix (without mls method)
    [Qsr,Qsc,Qsv,qc]    = assmbl_ful_q(iele,nnpe,np,qe,Qsr,Qsc,Qsv,qc);
    
    disp(['L2 fit + eq. constraint: ' num2str(iele/numel*100) '% is complete!'])
end

% free surface condition
% Calculate the constant matrices
bigNsurf    = bigNsurfmat(nnps,nqpts,ssfac);
% For each prescribed free surface
for   ieles = 1:1:numels
    % Surface normal and jacobian
    [n rjs] = surfjac(ieles,nnps,nqpts,dnda,dndb,nps,x,y,z);
    
    % Form the surface "stiffness matrix"
    % Assemble the free surface condition to the corresponding location
    % in the overall stiffness matrix
    sek = freesurf_residual(nnps,nqpts,swt,bigNsurf,n,rjs,beta);
    [Ksr,Ksc,Ksv,kc]                = assmblsurf_ful(ieles,nnps,nps,sek,Ksr,Ksc,Ksv,kc);
    
    disp(['Applying free surface constraint: ' num2str(ieles/numels*100) '% is complete!']);
end

% symmetry surface condition (in-plane surface tractions are free)
% For each prescribed free surface
kc_sym  = kc;
for   ieles = 1:1:numelsyms
    % Surface normal and jacobian
    [n rjs] = surfjac(ieles,nnps,nqpts,dnda,dndb,npsyms,x,y,z);
    
    % Form the surface "stiffness matrix"
    % Assemble the free surface condition to the corresponding location
    % in the overall stiffness matrix
    sek = symsurf_residual(nnps,nqpts,swt,bigNsurf,n,rjs,alpha);
    [Ksr,Ksc,Ksv,kc]                    = assmblsurf_ful(ieles,nnps,npsyms,sek,Ksr,Ksc,Ksv,kc);
    
    disp(['Applying symmetry surface constraint: ' num2str(ieles/numelsyms*100) '% is complete!'])
end

K   = sparse(Ksr,Ksc,Ksv);
Q   = sparse(Qsr,Qsc,Qsv);