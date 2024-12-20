clear all
close all
clc

% addpath(genpath('/home/beams/PARKJS/matlab/matlab_tools'));
addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER INPUT
wsname  = 'wscub3x';      % workspace name
% pname{1}    = '/net/wolf/data/tomo1/krause_mar23_midas/ff/LayerNr_1/';
pname{1}    = './analysis_cal_steel_ff_001052';
% pname{1}    = './analysis_Au_cubes_ff_001047_after_Tx_optimization';
% pname{1}    = './analysis_shrus_Al2O3_pyramid_planar_ff_001701_pass1';

fname   = 'Grains.csv';

% FILTERS - CHECK LINE 35 TO SEE WHICH ONE IS ON
Thresh_Completeness = 0.7;
Thresh_GrainRadius  = 50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Load workspace for fundamental region.
load(wsname);
eval(['ws = ', wsname, ';']);
clear(wsname)

% Load MIDAS results
for ii = 1:1:1
    pfname{ii}  = fullfile(pname{ii}, fname);
end

microstructure  = parseGrainData_MultiLayer_ff(pname, ws.frmesh.symmetries, ...
    'CrdSystem', 'APS', ...
    'LabToSample', 0, ...
    'C_xstal', BuildElasticityMatrix([176 224 97 61 51]*1000, 'Symmetry', 'hexagonal'), ...
    'OffsetValue', 100, ...
    'OutputReflectionTable', true);

Grains  = microstructure.grains;

numpts  = length(Grains);
wts     = ones(1, numpts);

% THRESHOLDING BY COMPLETENESS
idx_Completeness    = [Grains.Completeness] >= Thresh_Completeness;
idx_MeanRadius      = [Grains.GrainRadius] >= Thresh_GrainRadius;

idx = find(idx_Completeness);

grainID = [Grains(idx).GrainID]';
xyz         = [Grains(idx).COM]';
xyz(:,2)	= -xyz(:,2);
rod     = [Grains(idx).rod];
cidx    = [Grains(idx).Completeness]';
quat    = [Grains(idx).quat];
GrainRad    = [Grains(idx).GrainRadius]';
lattprm     = [Grains(idx).lattprms];
vm          = [Grains(idx).StressFab_vm]';
h           = [Grains(idx).StressFab_h];
StressFab_dev   = [Grains(idx).StressFab_d];
sxx_dev     = StressFab_dev(1,:)';
syy_dev     = StressFab_dev(2,:)';
szz_dev     = StressFab_dev(3,:)';
sxy_dev     = StressFab_dev(4,:)'./sqrt(2);
sxz_dev     = StressFab_dev(5,:)'./sqrt(2);
syz_dev     = StressFab_dev(6,:)'./sqrt(2);
StressFab       =  [Grains(idx).StressFab];
sxx         = StressFab(1,:)';
syy         = StressFab(2,:)';
szz         = StressFab(3,:)';
sxy         = StressFab(4,:)'./sqrt(2);
sxz         = StressFab(5,:)'./sqrt(2);
syz         = StressFab(6,:)'./sqrt(2);

% hdr_data    = 'x,y,z,GrainRad,sxx,syy,szz,sxy,sxz,syz,hs,vms,ds,completeness\n';
% pfname  = 'test.csv';
% fid = fopen(pfname, 'w');
% fprintf(fid, hdr_data);
% fclose(fid);
% dlmwrite('test.csv', [xyz GrainRad sxx syy szz sxy sxz syz h vm d cidx], 'delimiter', ',', '-append')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOTS IN PHYSICAL SPACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PLOT COM / ONE COLOR
figure,
if strcmpi( Grains(1).CrdSys, 'ESRF')
    scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 30, 'filled', 'b')
    xlabel('z : +=along beam (um)'); ylabel('x : +=OB (um)'); zlabel('y : +=UP (um)')
elseif strcmpi( Grains(1).CrdSys, 'APS')
    scatter3(xyz(:,1), xyz(:,3), xyz(:,2), 30, 'filled', 'b')
    xlabel('x : +=OB (um)'); ylabel('z : +=along beam (um)'); zlabel('y : +=UP (um)')
end
grid on; axis equal
title('COM of found grains')

%%%% PLOT COM / COMPLETENESS AS COLOR
figure, 
if strcmpi( Grains(1).CrdSys, 'ESRF')
    scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 30, cidx, 'filled')
    xlabel('z : +=along beam (um)'); ylabel('x : +=OB (um)'); zlabel('y : +=UP (um)')
elseif strcmpi( Grains(1).CrdSys, 'APS')
    scatter3(xyz(:,1), xyz(:,3), xyz(:,2), 30, cidx, 'filled')
    xlabel('x : +=OB (um)'); ylabel('z : +=along beam (um)'); zlabel('y : +=UP (um)')
end
colorbar vert; caxis([0.5 1.0])
grid on; axis equal; colormap jet
title('COM of found grains // colors denote completeness')

%%%% PLOT COM / VM AS COLOR
figure,
if strcmpi( Grains(1).CrdSys, 'ESRF')
    scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 30, vm, 'filled')
    xlabel('z : +=along beam (um)'); ylabel('x : +=OB (um)'); zlabel('y : +=UP (um)')
elseif strcmpi( Grains(1).CrdSys, 'APS')
    scatter3(xyz(:,1), xyz(:,3), xyz(:,2), 30, vm, 'filled')
    xlabel('x : +=OB (um)'); ylabel('z : +=along beam (um)'); zlabel('y : +=UP (um)')
end
grid on; axis equal; colormap jet
colorbar vert; caxis([min(vm) max(vm)])
title('COM of found grains // colors denote VM')

%%%% PLOT COM / H AS COLOR
figure,
if strcmpi( Grains(1).CrdSys, 'ESRF')
    scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 30, h, 'filled')
    xlabel('z : +=along beam (um)'); ylabel('x : +=OB (um)'); zlabel('y : +=UP (um)')
elseif strcmpi( Grains(1).CrdSys, 'APS')
    scatter3(xyz(:,1), xyz(:,3), xyz(:,2), 30, h, 'filled')
    xlabel('x : +=OB (um)'); ylabel('z : +=along beam (um)'); zlabel('y : +=UP (um)')
end
grid on; axis equal; colormap jet
colorbar vert; caxis([min(h) max(h)])
title('COM of found grains // colors denote H')

titlestr{1} = 'COM of found grains // colors denote sxx';
titlestr{2} = 'COM of found grains // colors denote syy';
titlestr{3} = 'COM of found grains // colors denote szz';
titlestr{4} = 'COM of found grains // colors denote sxy';
titlestr{5} = 'COM of found grains // colors denote sxz';
titlestr{6} = 'COM of found grains // colors denote syz';
for i = 1:1:6
    if i < 4
        sij = StressFab(i,:);
    else
        sij = StressFab(i,:)./sqrt(2);
    end
    %%%% PLOT COM / sij AS COLOR
    figure, 
    if strcmpi( Grains(1).CrdSys, 'ESRF')
        scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 30, sij, 'filled')
        xlabel('z : +=along beam (um)'); ylabel('x : +=OB (um)'); zlabel('y : +=UP (um)')
    elseif strcmpi( Grains(1).CrdSys, 'APS')
        scatter3(xyz(:,1), xyz(:,3), xyz(:,2), 30, sij, 'filled')
        xlabel('x : +=OB (um)'); ylabel('z : +=along beam (um)'); zlabel('y : +=UP (um)')
    end
    grid on; axis equal; colormap jet
    title(titlestr{i})
    colorbar vert; caxis([min(sij) max(sij)])
end