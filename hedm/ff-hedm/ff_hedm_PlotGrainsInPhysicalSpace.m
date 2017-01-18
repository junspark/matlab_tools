clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER INPUT
wsname  = 'wscub4x';      % workspace name
pname   = './examples/';
fname   = 'Grains.csv';

% ROTATION MATRIX TAKING VECTOR IN LAB FRAME TO SAMPLE FRAME
% NECESSARY TO GET THE ORIENTATION OF CRYSTALS WITH RESPECT TO SAMPLE FRAME
RLab2Sam    = eye(3,3);

% FILTERS - CHECK LINE 35 TO SEE WHICH ONE IS ON
Thresh_Completeness = 0.7;
Thresh_GrainRadius  = 50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Load workspace for fundamental region.
load(wsname);
eval(['ws = ', wsname, ';']);
clear(wsname)

% Load MIDAS results
pfname  = fullfile(pname, fname);
Grains  = parseGrainData(pfname, ws.frmesh.symmetries, ...
    'CrdSystem', 'ESRF', ...
    'LabToSample', 0, ...
    'C_xstal', BuildElasticityMatrix([176 124 81.7], 'Symmetry', 'cubic'));
numpts  = length(Grains);
wts     = ones(1, numpts);

% THRESHOLDING BY COMPLETENESS
idx_Completeness    = [Grains.Completeness] >= Thresh_Completeness;
idx_MeanRadius      = [Grains.GrainRadius] >= Thresh_GrainRadius;

idx = find(idx_Completeness);

grainID = [Grains(idx).GrainID]';
xyz     = [Grains(idx).COM]';
rod     = [Grains(idx).rod];
cidx    = [Grains(idx).Completeness];
quat    = [Grains(idx).quat];
GrainRad    = [Grains(idx).GrainRadius];
lattprm     = [Grains(idx).lattprms];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOTS IN PHYSICAL SPACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PLOT COM / ONE COLOR
figure, scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 30, 'filled', 'b')
grid on; axis square
xlabel('z : +=along beam (um)'); ylabel('x : +=OB (um)'); zlabel('y : +=UP (um)')
title('COM of found grains')

%%%% PLOT COM / COMPLETENESS AS COLOR
figure, scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 30, cidx, 'filled')
grid on; axis square
colorbar vert; caxis([0.5 1])
xlabel('z : +=along beam (um)'); ylabel('x : +=OB (um)'); zlabel('y : +=UP (um)')
title('COM of found grains // colors denote completeness')