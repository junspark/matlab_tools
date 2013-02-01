%%%% THIS DOES NOT WORK AT THE MOMENT BECAUSE WE ARE USING MESHLESS METHOD
clear all
close all
clc

load inputdata.A1.mat
load Residual_Stress_3D_MLSEL_v1/meshdata.mat

PFNAME_Re   = fullfile(inputdata.PNAME_SAM, inputdata.FNAME_Re);
PFNAME_RANK = fullfile(inputdata.PNAME_SAM, inputdata.FNAME_RANK);

load(PFNAME_Re)     % Re
load(PFNAME_RANK)   % Rk

nx  = length(inputdata.x_grid);
ny  = length(inputdata.y_grid);
nz  = length(inputdata.z_grid);
npf = size(inputdata.hkls,1);

% FIND MAX AND MIN OF Re
max_Re  = 0;
min_Re  = 0;
for i = 1:1:npf
    if max(Re{i}(:,4)) > max_Re
        max_Re  = max(Re{i}(:,4));
    end
    if min(Re{i}(:,4)) < min_Re
        min_Re  = min(Re{i}(:,4));
    end
    
end

%%% PLOT Re FOR EACH SPF
for i = 1:1:npf
    f   = figure(i);
    c   = Re{i}(:,4);
    return
    PlotDV(f, np, x, y, z, c);
    colorbar('location', 'SouthOutside')
    caxis([min_Re max_Re])
    
    titlestr    = ['R_{e} - \{', ...
        num2str(inputdata.hkls(i,1)), ...
        num2str(inputdata.hkls(i,2)), ...
        num2str(inputdata.hkls(i,3)), ...
        '\}'];
    title(titlestr)
    view([45 26])
    axis equal tight   
end

f   = figure(npf+1);
c   = Re{size(npf,1)+1}(:,4);
PlotDV(f, np, x, y, z, c);
colorbar('location', 'SouthOutside')
caxis([min_Re max_Re])
title('R_{e} - ALL SPF')
view([45 26])
axis equal tight

f   = figure(npf+2);
c   = Rk(:,4);
PlotDV(f, np, x, y, z, c, 'ShowMesh', 'on');
colorbar('location', 'SouthOutside')
title('RANK - INV MATRIX')
view([45 26])
axis equal tight