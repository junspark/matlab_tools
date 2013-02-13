clear all
close all
clc

%%% INPUT DATA GENERATED IN STEP0
InputData   = load('example_rsinput.mat');

%%% DV INFO
PFNAME_GRID = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_GRID);
GridData    = load(PFNAME_GRID);

%%% MESH INFO
PFNAME_MESH = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_MESH);
MeshData    = load(PFNAME_MESH);

PFNAME_Re   = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_Re);
load(PFNAME_Re)

%%% 
% MACRO STRESS
np  = MeshData.np;
x   = MeshData.x;
y   = MeshData.y;
z   = MeshData.z;

numnp   = MeshData.numnp;
numel   = size(np,1);
ncmap   = 256;
cmap    = jet(ncmap);
for i = 1:1:InputData.nhkls
    % define colors
    Data    = Re{i}(:,4);
    DataRange(1)    = min(Data);
    DataRange(2)    = max(Data);
    
    dData   = DataRange(2) - DataRange(1);
    DVColor = round((ncmap-1)*(Data - DataRange(1))./dData) + 1;
    
    %%% PLOT MESH FIRST
    figure(i)
    hold on
    
    for j = 1:1:numel
        n1  = np(j,1);
        n2  = np(j,2);
        n3  = np(j,3);
        n4  = np(j,4);
        n5  = np(j,5);
        n6  = np(j,6);
        n7  = np(j,7);
        n8  = np(j,8);
        
        xx  = [x(n1) x(n2) x(n3) x(n4)...
            x(n5) x(n6) x(n7) x(n8)];
        yy  = [y(n1) y(n2) y(n3) y(n4)...
            y(n5) y(n6) y(n7) y(n8)];
        zz  = [z(n1) z(n2) z(n3) z(n4)...
            z(n5) z(n6) z(n7) z(n8)];
        
        line([xx(1:4) xx(1)],[yy(1:4) yy(1)],[zz(1:4) zz(1)],'Color','k');
        line([xx(5:8) xx(5)],[yy(5:8) yy(5)],[zz(5:8) zz(5)],'Color','k');
        line([xx(1) xx(5)],[yy(1) yy(5)],[zz(1) zz(5)],'Color','k');
        line([xx(2) xx(6)],[yy(2) yy(6)],[zz(2) zz(6)],'Color','k');
        line([xx(3) xx(7)],[yy(3) yy(7)],[zz(3) zz(7)],'Color','k');
        line([xx(4) xx(8)],[yy(4) yy(8)],[zz(4) zz(8)],'Color','k');
        plot3(xx, yy, zz, 'k.')
        hold on
    end
    
    for j = 1:1:GridData.ndv
        if DVColor(j) > 0 & DVColor(j) <= ncmap
            plot3(...
                Re{i}(j,1), Re{i}(j,2), Re{i}(j,3), ...
                'Marker', 'o', ...
                'MarkerFaceColor', cmap(DVColor(j), :), ...
                'MarkerEdgeColor', cmap(DVColor(j), :), ...
                'MarkerSize', 10);
            hold on
        end
    end
    view([45 30])
    axis equal tight off
    
    cb  = colorbar;
    set(cb, 'FontSize', 16, 'FontWeight', 'bold', ...
        'Location', 'EastOutside')
    caxis(DataRange);
    titlestr    = ['{hk.l} = {', ...
        num2str(InputData.hkls(i,1)), ...
        num2str(InputData.hkls(i,2)), ...
        num2str(InputData.hkls(i,3)), ...
        '}'];
    title(titlestr)
end