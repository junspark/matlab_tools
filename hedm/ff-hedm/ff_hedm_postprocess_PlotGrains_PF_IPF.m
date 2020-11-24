clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INSTALLS MTEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run('/home/beams/S1IDUSER/mnt/s1b/__eval/mtex-4.4.0/startup_mtex.m')

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% workspace name
wsname  = 'wscub4x';

%%% MIDAS RESULT PATH
pname_root  = '/home/beams/S1IDUSER/mnt/orthros/internal_aug18_midas/ff/hedm_resolution';
pname       = fullfile(pname_root, 'ss_sam_ff3_Layer8_Analysis_Time_2018_09_11_12_27_14');

%%% IDENTIFY GRAIN OF INTEREST USING GRAINID
gid         = 5414;

%%% MTEX RELATED PARAMETERS
%%% CRYSTAL SYMMETRY IN MTEX CONVENTION
cs  = crystalSymmetry('Td', [3.662 3.662 3.662], 'mineral', 'Iron', 'color', 'light blue');

% crystal symmetry
CS = {... 
    'notIndexed', ...
    cs};

%%% CRYSTAL SYMMETRY IN MTEX CONVENTION
ss  = specimenSymmetry('orthorhombic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOAD DATA AND WORKSPACES
grain_map   = parseGrainData_OneLayer(pname, CubSymmetries, ...
    'CrdSystem', 'APS', ...
    'LabToSample', 0, ...
    'C_xstal', nan, ...
    'ComputeSelfMisoTable', false, ...
    'ComputeSelfDistTable', false, ...
    'OutputReflectionTable', false, ...
    'NumFrames', 1440, ...
    'Technique', 'ff-midas');

%%% Import Grains.csv file for EBSD Data
pfname  = fullfile(pname, 'Grains.csv');

%%% Import the Data
% create an EBSD variable containing the data
ebsd = loadEBSD(pfname, CS, 'interface', 'generic', ...
    'ColumnNames', { 'x' 'y' 'z' 'Weight' 'ConfidenceIndex' 'Phase' 'phi1' 'Phi' 'phi2'}, ...
    'Columns', [11 12 13 23 24 44 45 46 47], 'Bunge');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% POLE FIGURE PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
plotPDF(ebsd('Iron').orientations,Miller(1,0,0,ebsd('Iron').CS))

%%% UNINSTALL MTEX IN CASE OF CONFLICTS
% run('/home/beams/S1IDUSER/mnt/s1b/__eval/mtex-4.4.0/uninstall_mtex.m')

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
