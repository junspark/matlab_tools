% This code is written to generate the mesh and the boundary conditions for
% the shrink-fit-disk problem
clear
close all
clc
%
% Inputs (diffraction data centers)
% inputs_KM
% inputs
inputs_A1234
load('/home/jspark/Documents/work/prjResidualStress/3D.MultiscaleCode2/inputdata.A1.mat')
PFNAME_GRID = fullfile(inputdata.PNAME_SAM, inputdata.FNAME_GRID);
load(PFNAME_GRID)

% inputs
%
% Calculate DV stresses using the analytical solution
% [ndv,dvc_x,dvc_y,dvc_z]=findnodecoordinates(R,t,alpha);
%
% Values of R, t, and alpha at the nodes'
alpha   = [90:10:180];
dalpha  = 0;
t       = -0.6:0.2:0.6;
dt      = 0.05;
[Rnodes,tnodes,alphanodes]  = findnodalvalues(R,t,dt,alpha,dalpha);

% Find the coordinates of the nodes
[numnp,x,y,z]   = findnodecoordinates(Rnodes,tnodes,alphanodes);

for i = 1:1:numnp
   plot3(x(i), y(i), z(i), 'k.') 
   hold on
end
% Find the connectivity matrix
[np]            = findconnectivity(Rnodes,tnodes,alphanodes);

numel=size(np,1);
for i = 1:1:numel
    for j = 1:1:8
        plot3(x(np(i,j)), y(np(i,j)), z(np(i,j)), 'r.')
    end
end

%
% Apply boundary conditions
% traction free : nps numels
% shear free : npsyms, numelsyms
[nps,numels,npsyms,numelsyms]=applyBCs(Rnodes,tnodes,alphanodes);
%
% Plot the mesh
plotmesh
%
%
% Mesh properties (8-node brick)
nnps    = 4;
meltyp  = 8;
numeq   = 3*numnp;
nnpe    = meltyp;

% Save the meshdata
save meshdata numnp numeq numel numels numelsyms meltyp nnpe nnps np nps npsyms x y z

% for i = 1:1:numel
%     np(i,:);
%     plot3(x(np(i,1)), y(np(i,1)), z(np(i,1)), 'r.')
%     plot3(x(np(i,2)), y(np(i,2)), z(np(i,2)), 'g.')
%     plot3(x(np(i,3)), y(np(i,3)), z(np(i,3)), 'b.')
%     plot3(x(np(i,4)), y(np(i,4)), z(np(i,4)), 'm.')
%     plot3(x(np(i,5)), y(np(i,5)), z(np(i,5)), 'rs')
%     plot3(x(np(i,6)), y(np(i,6)), z(np(i,6)), 'gs')
%     plot3(x(np(i,7)), y(np(i,7)), z(np(i,7)), 'bs')
%     plot3(x(np(i,8)), y(np(i,8)), z(np(i,8)), 'ms')
%     plot3(mean(x(np(i,:))), mean(y(np(i,:))), mean(z(np(i,:))), 'y.')
% %     pause
% end

for i = 1:1:numelsyms
    npsyms(i,:);
    for j = 1:1:4
        plot3(x(npsyms(i,j)), y(npsyms(i,j)), z(npsyms(i,j)), 'ro')
        % pause
    end
end

for i = 1:1:numels
    nps(i,:);
    for j = 1:1:4
        plot3(x(nps(i,j)), y(nps(i,j)), z(nps(i,j)), 'gs')
        % pause
    end
end

for i = 1:1:ndv
    plot3(dvc_x(i), dvc_y(i), dvc_z(i), 'b^')
    % pause(0.5)
end
disp('done')
