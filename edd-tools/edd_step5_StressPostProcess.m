clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSTRUMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NUMBER OF POINTS IN X,Y,Z (GET FROM PYTHON FILE)
NX  = 15; NY  = 1; NZ  = 7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATERIAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E   = 200000;       % [MPa]
nu  = 0.3;          % [-]
pkid_ave    = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PAR FILE DESIGNATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pname_pypar     = './strain-examples/';
fname_pypar     = 'mach_feb15_TOA_7_Lap7.pypar';
pfname_pypar    = fullfile(pname_pypar, fname_pypar);
pardata         = ReadPythonParFile(pfname_pypar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATH WHERE DATA FILES LIVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pname_data  = './strain-examples/Lap7/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILE NAME FOR STRAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname_strain_v  = 'mach_feb15_strain_v.csv';
fname_strain_h  = 'mach_feb15_strain_h.csv';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XS  = pardata.samX;
YS  = pardata.samY;
ZS  = pardata.samZ;

xs  = reshape(XS', NX, NY, NZ);
ys  = reshape(YS', NX, NY, NZ);
zs  = reshape(ZS', NX, NY, NZ);

numDV   = length(pardata.day);

strain_v    = csvread(fname_strain_v);
strain_h    = csvread(fname_strain_h);
    
epsilon_vv  = mean(strain_v(:,pkid_ave + 3), 2);
epsilon_hh  = mean(strain_h(:,pkid_ave + 3), 2);

sigma_vv    = E/(1 + nu)/(1 - 2*nu) * ((1 - nu) * epsilon_vv + nu * epsilon_hh);
sigma_hh    = E/(1 + nu)/(1 - 2*nu) * (nu * epsilon_vv + (1 - nu) * epsilon_hh);
sigma_oo    = E/(1 + nu)/(1 - 2*nu) * (nu * epsilon_vv + nu * epsilon_hh);

figure(1)
% set(gcf, 'Position', [120 560 1640 420])
subplot(1,2,1)
scatter3(XS, YS, ZS, 50, epsilon_vv, 'filled')
colorbar vert
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
caxis([-1e-3 2e-3])
title('epsilon_{vv}')

subplot(1,2,2)
scatter3(XS, YS, ZS, 50, epsilon_hh, 'filled')
colorbar vert
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
caxis([-1e-3 2e-3])
title('epsilon_{hh}')

figure(2)
% set(gcf, 'Position', [120 560 1640 420])
subplot(1,3,1)
scatter3(XS, YS, ZS, 50, sigma_vv, 'filled')
colorbar vert
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
caxis([-300 300])
title('sigma_{vv} (assuming plane strain / no shear strain)')

subplot(1,3,2)
scatter3(XS, YS, ZS, 50, sigma_hh, 'filled')
colorbar vert
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
caxis([-300 300])
title('sigma_{hh} (assuming plane strain / no shear strain)')

subplot(1,3,3)
scatter3(XS, YS, ZS, 50, sigma_oo, 'filled')
colorbar vert
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
caxis([-300 300])
title('sigma_{oo} (assuming plane strain / no shear strain)')