clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INSTALLS MTEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run('/home/beams/S1IDUSER/mnt/s1b/__eval/mtex-4.4.0/startup_mtex.m')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% workspace name
wsname_cub  = 'wscub4x';

%%% MIDAS RESULT PATH
pname_root  = '/home/beams/S1IDUSER/mnt/orthros/internal_aug18_midas/ff/hedm_resolution';
pname       = fullfile(pname_root, 'ss_sam_ff3_Layer8_Analysis_Time_2018_09_11_12_27_14');
gid         = 5414;

%%% CRYSTAL SYMMETRY IN MTEX CONVENTION
cs  = crystalSymmetry('Td');
CS  = {...
    crystalSymmetry('Td', [3.6 3.6 3.6], 'mineral', 'Fe', 'color', 'light blue'), ...
    };

%%% CRYSTAL SYMMETRY IN MTEX CONVENTION
ss  = specimenSymmetry('orthorhombic');
SS  = specimenSymmetry('orthorhombic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(wsname_cub);
eval(['wscub    = ', wsname_cub, ';']);
clear(wsname_cub)

grain_map   = parseGrainData_OneLayer(pname, CubSymmetries, ...
    'CrdSystem', 'APS', ...
    'LabToSample', 0, ...
    'C_xstal', nan, ...
    'ComputeSelfMisoTable', 1, ...
    'ComputeSelfDistTable', 1, ...
    'OutputReflectionTable', false, ...
    'NumFrames', 1440, ...
    'Technique', 'ff-midas');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GRAIN OF INTEREST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GrainID = [grain_map.grains.GrainID];
idx     = find(GrainID == gid);

rod     = grain_map.grains(idx).rod;
bunge   = BungeOfRMat(grain_map.grains(idx).R, 'degrees');
ori     = orientation('Euler', bunge(1)*degree, bunge(2)*degree, bunge(3)*degree, cs, ss);

%%% ORI PLOTS
figure, 
PlotFRPerimeter('cubic');
axis square off
hold on
scatter3(rod(1), rod(2), rod(3), 'o', 'r', 'filled')
hold off

%%% IPF PLOTS
figure,
plotIPDF(ori, yvector, cs, 'antipodal')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ALL GRAINS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
quat    = [grain_map.grains(:).quat];
rod     = [grain_map.grains(:).rod];
RMat    = RMatOfQuat(quat);
bunge   = BungeOfRMat(RMat, 'degrees');
ori     = orientation('Euler', bunge(1,:)*degree, bunge(2,:)*degree, bunge(3,:)*degree, cs, ss);

%%% ORI PLOTS
figure, 
PlotFRPerimeter('cubic');
axis square off
hold on
scatter3(rod(1,:), rod(2,:), rod(3,:), 'o', 'r', 'filled')
hold off

%%% IPF PLOTS
figure,
plotIPDF(ori, yvector, cs, 'antipodal')
hold off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% POLE FIGURE PLOTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Import Script for EBSD Data
%
% This script was automatically created by the import wizard. You should
% run the whoole script or parts of it in order to import your data. There
% is no problem in making any changes to this script.

%%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [3.662 3.662 3.662], 'mineral', 'Iron', 'color', 'light blue')};

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%%% Specify File Names

% path to files
pname = '/home/s1b/__eval/projects_parkjs/ff_plotter';

% which files to be imported
fname = [pname '/Grains.csv'];

%%% Import the Data

% create an EBSD variable containing the data
ebsd = loadEBSD(fname,CS,'interface','generic',...
  'ColumnNames', { 'x' 'y' 'z' 'Weight' 'ConfidenceIndex' 'Phase' 'phi1' 'Phi' 'phi2'}, 'Columns', [11 12 13 23 24 44 45 46 47], 'Bunge');
 
plotPDF(ebsd('Iron').orientations,Miller(1,0,0,ebsd('Iron').CS))

%%% UNINSTALLS MTEX IN CASE OF CONFLICTS
run('/home/beams/S1IDUSER/mnt/s1b/__eval/mtex-4.4.0/uninstall_mtex.m')

% c   = [1 1 1]';
% cu  = UnitVector(c);
% for i = 1:1:size(ws.frmesh.symmetries,2)
%     R   = RMatOfQuat(ws.frmesh.symmetries(:,i));
%     ci(:,i) = R*cu;
% end
% 
% s   = [];
% for i = 1:1:size(quat,2)
%     R   = RMatOfQuat(quat(:,i));
%     s   = [s R*ci];
% end
% 
% PlotSPF(s', ones(size(s,2),1), 'Title', sprintf('{ %d%d%d } in sample frame', c(1), c(2), c(3)))
