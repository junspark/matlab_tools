clear all
close all
clc

%%% INPUT DATA GENERATED IN STEP0
InputData   = load('example_rsinput.mat');

%%% FILE NAME FOR SOLUTION
FNAME_SOL   = 'Solution_l1_1_l2_0.1_l3_0.1_k_1.5.mat';
PFNAME_SOL  = fullfile(InputData.PNAME_SOLUTION, FNAME_SOL);
SOL = load(PFNAME_SOL);
SOL = SOL.SOL2;

%%% MESH INFO
PFNAME_MESH = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_MESH);
MeshData    = load(PFNAME_MESH);

%%% DV INFO
PFNAME_GRID = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_GRID);
GridData    = load(PFNAME_GRID);

%%% MICRO MATRIX
PFNAME_A    = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_A);
load(PFNAME_A)
PFNAME_N    = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_N);
load(PFNAME_N)

%%% MACRO MATRICES
PFNAME_KQ    = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_KQ);
load(PFNAME_KQ)

%%% 
% MACRO STRESS
numnp       = MeshData.numnp;
Stress_h    = SOL(1:numnp*6,1).*1000000;  %%% TPa <> MPa

%%%
% MACRO STRESS
Stress_h11  = Stress_h(1:6:numnp*6-5);
Stress_h22  = Stress_h(2:6:numnp*6-4);
Stress_h33  = Stress_h(3:6:numnp*6-3);
Stress_h23  = Stress_h(4:6:numnp*6-2);
Stress_h13  = Stress_h(5:6:numnp*6-1);
Stress_h12  = Stress_h(6:6:numnp*6-0);

Stress_hVM  = sqrt((...
    (Stress_h11 - Stress_h22).^2 + ...
    (Stress_h22 - Stress_h33).^2 + ...
    (Stress_h33 - Stress_h11).^2 + ...
    6*Stress_h12.*Stress_h12 + ...
    6*Stress_h13.*Stress_h13 + ...
    6*Stress_h23.*Stress_h23)./2);

f   = figure(1000);
c   = zeros(numnp,1);
hold on
PlotRSMesh(f, MeshData.np, MeshData.x, MeshData.y, MeshData.z, c, ...
    'ShowMesh', 'on', ...
    'DataRange', [-1 +1])

f   = figure(1);
subplot(2,3,1)
hold on
PlotRSMesh(f, MeshData.np, MeshData.x, MeshData.y, MeshData.z, Stress_h11, ...
    'ShowMesh', 'on', ...
    'DataRange', [-100 100], ...
    'Title', '\Sigma_{11}', ...
    'Colorbar', 'on')

f   = figure(1);
subplot(2,3,2)
hold on
PlotRSMesh(f, MeshData.np, MeshData.x, MeshData.y, MeshData.z, Stress_h22, ...
    'ShowMesh', 'on', ...
    'DataRange', [-100 100], ...
    'Title', '\Sigma_{22}', ...
    'Colorbar', 'on')

f   = figure(1);
subplot(2,3,3)
hold on
PlotRSMesh(f, MeshData.np, MeshData.x, MeshData.y, MeshData.z, Stress_h33, ...
    'ShowMesh', 'on', ...
    'DataRange', [-50 900], ...
    'Title', '\Sigma_{22}', ...
    'Colorbar', 'on')

f   = figure(1);
subplot(2,3,4)
hold on
PlotRSMesh(f, MeshData.np, MeshData.x, MeshData.y, MeshData.z, Stress_h23, ...
    'ShowMesh', 'on', ...
    'DataRange', [-100 100], ...
    'Title', '\Sigma_{23}', ...
    'Colorbar', 'on')

f   = figure(1);
subplot(2,3,5)
hold on
PlotRSMesh(f, MeshData.np, MeshData.x, MeshData.y, MeshData.z, Stress_h13, ...
    'ShowMesh', 'on', ...
    'DataRange', [-100 100], ...
    'Title', '\Sigma_{13}', ...
    'Colorbar', 'on')

f   = figure(1);
subplot(2,3,6)
hold on
PlotRSMesh(f, MeshData.np, MeshData.x, MeshData.y, MeshData.z, Stress_h12, ...
    'ShowMesh', 'on', ...
    'DataRange', [-100 100], ...
    'Title', '\Sigma_{12}', ...
    'Colorbar', 'on')

f   = figure(2);
hold on
PlotRSMesh(f, MeshData.np, MeshData.x, MeshData.y, MeshData.z, Stress_hVM, ...
    'ShowMesh', 'on', ...
    'DataRange', [-50 900], ...
    'Title', '\Sigma_{VM}', ...
    'Colorbar', 'on')