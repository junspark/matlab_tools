clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER INPUT
wsname  = 'wshex3x';      % workspace name
pname{1}    = '/net/wolf/data/tomo1/broderick_dec16_MIDAS/ff/HeatHTNS9/S8/HeatHTNS9_S8_FF1_Layer1_Analysis_Time_2017_01_19_17_08_22';
pname{2}    = '/net/wolf/data/tomo1/broderick_dec16_MIDAS/ff/HeatHTNS9/S8/HeatHTNS9_S8_FF1_Layer2_Analysis_Time_2017_01_19_18_34_06';
pname{3}    = '/net/wolf/data/tomo1/broderick_dec16_MIDAS/ff/HeatHTNS9/S8/HeatHTNS9_S8_FF1_Layer3_Analysis_Time_2017_01_19_19_43_57';
pname{4}    = '/net/wolf/data/tomo1/broderick_dec16_MIDAS/ff/HeatHTNS9/S8/HeatHTNS9_S8_FF1_Layer4_Analysis_Time_2017_01_19_21_49_05';
pname{5}    = '/net/wolf/data/tomo1/broderick_dec16_MIDAS/ff/HeatHTNS9/S8/HeatHTNS9_S8_FF1_Layer5_Analysis_Time_2017_01_20_00_22_23';
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
for ii = 1:1:5
    pfname{ii}  = fullfile(pname{ii}, fname);
end

Grains	= parseGrainData(pfname, ws.frmesh.symmetries, ...
    'CrdSystem', 'APS', ...
    'LabToSample', 0, ...
    'C_xstal', BuildElasticityMatrix([176 224 97 61 51]*1000, 'Symmetry', 'hexagonal'), ...
    'OffsetDirection', 'y', ...
    'OffsetValue', 100);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STATISTICAL PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HISTOGRAM OF GRAIN SIZES NORMALIZED BY MAX GRAIN SIZE
figure, 
hist(GrainRad./max(GrainRad), 20)
xlabel('relative grain radius (-)')
ylabel('number of grains (-)')
title(sprintf('Max grain size : %5.0f (micron)', max(GrainRad)))

figure, 
subplot(2,3,1)
hist(lattprm(1,:))
xlabel('a (Angstrom)')
ylabel('number of grains (-)')
% title(sprintf('a0 = %5.4f A', a0))
view([0 90])
% axis([3.58 3.61 0 150])
grid on

subplot(2,3,2)
hist(lattprm(2,:))
xlabel('b (Angstrom)')
ylabel('number of grains (-)')
% title(sprintf('a0 = %5.4f A', a0))
view([0 90])
% axis([3.58 3.61 0 150])
grid on

subplot(2,3,3)
hist(lattprm(3,:))
xlabel('c (Angstrom)')
ylabel('number of grains (-)')
% title(sprintf('a0 = %5.4f A', a0))
view([0 90])
% axis([3.58 3.61 0 150])
grid on

subplot(2,3,4)
hist(lattprm(4,:))
xlabel('\alpha (degrees)')
ylabel('number of grains (-)')
view([0 90])
% axis([89.7 90.3 0 150])
grid on

subplot(2,3,5)
hist(lattprm(5,:))
xlabel('\beta (degrees)')
ylabel('number of grains (-)')
view([0 90])
% axis([89.7 90.3 0 150])
grid on

subplot(2,3,6)
hist(lattprm(6,:))
xlabel('\gamma (degrees)')
ylabel('number of grains (-)')
view([0 90])
% axis([89.7 90.3 0 150])
grid on