clear all
close all
clc

%%% INPUT DATA GENERATED IN STEP0
InputData   = load('example_rsinput.mat');

dvc_x   = InputData.dvx;
dvc_y   = InputData.dvy;
dvc_z   = InputData.dvz;

ndv = InputData.ndv;

PFNAME_GRID = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_GRID);

disp('generating DV grid file')
disp('file name is ...')
disp(InputData.FNAME_GRID)

plot3(dvc_x, dvc_y, dvc_z, '.')
axis square
title('diffraction volume locations')
save(PFNAME_GRID, 'dvc_x', 'dvc_y', 'dvc_z', 'ndv')