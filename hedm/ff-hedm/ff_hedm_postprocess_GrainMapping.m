clear all
close all
clc

%%%
wsname_cub  = 'wscub4x';      % workspace name
load(wsname_cub);
eval(['wscub    = ', wsname_cub, ';']);
clear(wsname_cub)
%%%

ang_res = 0.25; % deg
pos_res = 100;  % um
vol_diff_thresh = 0.02; %%% THIS IS RADIUS FILTER

pname_root  = '/home/beams/S1IDUSER/mnt/orthros/internal_aug18_midas/ff/hedm_resolution';
pname0      = fullfile(pname_root, 'ss_sam_ff3_Layer8_Analysis_Time_2018_09_11_12_27_14');
pname1      = fullfile(pname_root, 'ss_sam_ff3_Layer1_Analysis_Time_2018_09_11_10_37_29');

save_mapping_table  = 1;
pfname_mapping_table    = 'ouput.mat';

S0_grains   = parseGrainData_OneLayer(pname0, CubSymmetries, ...
    'CrdSystem', 'APS', ...
    'LabToSample', 0, ...
    'C_xstal', nan, ...
    'ComputeSelfMisoTable', 1, ...
    'ComputeSelfDistTable', 1, ...
    'OutputReflectionTable', false, ...
    'NumFrames', 1440, ...
    'Technique', 'ff-midas');

S1_grains   = parseGrainData_OneLayer(pname1, CubSymmetries, ...
    'CrdSystem', 'APS', ...
    'LabToSample', 0, ...
    'C_xstal', nan, ...
    'ComputeSelfMisoTable', 1, ...
    'ComputeSelfDistTable', 1, ...
    'OutputReflectionTable', false, ...
    'NumFrames', 1440, ...
    'Technique', 'ff-midas');

mapping_table   = hedm_TrackGrains(S0_grains.grains, S1_grains.grains, ...
    'ang_res', ang_res, 'pos_res', pos_res, ...
    'save_mapping_table', save_mapping_table, ...
    'mapping_table_fname', pfname_mapping_table);

mapping_table_data  = load(pfname_mapping_table);

vol_diff_pct    = abs(mapping_table_data.mapping_table(:,22) - mapping_table_data.mapping_table(:,23))./mapping_table_data.mapping_table(:,22);
idx_vol_diff    = vol_diff_pct < vol_diff_thresh;

numgrains               = size(mapping_table_data.dist_table,2);
num_vol_matching_grain	= sum(idx_vol_diff);

%%%
figure(1)
subplot(1,2,1)
scatter3(mapping_table_data.mapping_table(:,8), mapping_table_data.mapping_table(:,9), mapping_table_data.mapping_table(:,10), ...
    20, mapping_table_data.mapping_table(:,6), 'filled')
colormap jet
caxis([median(mapping_table_data.mapping_table(:,6), 'omitnan')-0.1 median(mapping_table_data.mapping_table(:,6), 'omitnan')+0.1])
view([0 0])
axis square tight
hcb = colorbar;
set(hcb, 'Location', 'southoutside')
hcb.Label.String    = 'misorientation (deg)';
title('without volume filter')
xlabel('X (\mum)')
zlabel('Z (\mum)')

subplot(1,2,2)
scatter3(mapping_table_data.mapping_table(idx_vol_diff,8), mapping_table_data.mapping_table(idx_vol_diff,9), mapping_table_data.mapping_table(idx_vol_diff,10), ...
    20, mapping_table_data.mapping_table(idx_vol_diff,6), 'filled')
colormap jet
caxis([median(mapping_table_data.mapping_table(:,6), 'omitnan')-0.1 median(mapping_table_data.mapping_table(:,6), 'omitnan')+0.1])
view([0 0])
axis square tight
hcb = colorbar;
set(hcb, 'Location', 'southoutside')
hcb.Label.String    = 'misorientation (deg)';
title('with volume filter')
xlabel('X (\mum)')
zlabel('Z (\mum)')

%%%
figure(10)
subplot(1,2,1)
scatter3(mapping_table_data.mapping_table(:,8), mapping_table_data.mapping_table(:,9), mapping_table_data.mapping_table(:,10), ...
    20, mapping_table_data.mapping_table(:,22), 'filled')
colormap jet
caxis([min(mapping_table_data.mapping_table(:,22))-0.1 max(mapping_table_data.mapping_table(:,22))+0.1])
view([0 0])
axis square tight
hcb = colorbar;
set(hcb, 'Location', 'southoutside')
hcb.Label.String    = 'grain radius (\mum)';
title('without volume filter')
xlabel('X (\mum)')
zlabel('Z (\mum)')

subplot(1,2,2)
scatter3(mapping_table_data.mapping_table(idx_vol_diff,8), mapping_table_data.mapping_table(idx_vol_diff,9), mapping_table_data.mapping_table(idx_vol_diff,10), ...
    20, mapping_table_data.mapping_table(idx_vol_diff,22), 'filled')
colormap jet
caxis([min(mapping_table_data.mapping_table(:,22))-0.1 max(mapping_table_data.mapping_table(:,22))+0.1])
view([0 0])
axis square tight
hcb = colorbar;
set(hcb, 'Location', 'southoutside')
hcb.Label.String    = 'grain radius (\mum)';
title('with volume filter')
xlabel('X (\mum)')
zlabel('Z (\mum)')

figure(11)
scatter3(mapping_table_data.mapping_table(~idx_vol_diff,8), mapping_table_data.mapping_table(~idx_vol_diff,9), mapping_table_data.mapping_table(~idx_vol_diff,10), ...
    20, mapping_table_data.mapping_table(~idx_vol_diff,22), 'filled')
colormap jet
caxis([min(mapping_table_data.mapping_table(:,22))-0.1 max(mapping_table_data.mapping_table(:,22))+0.1])
view([0 0])
axis square tight
title('COMs of missing grains')
hcb = colorbar;
set(hcb, 'Location', 'eastoutside')
hcb.Label.String    = 'grain radius (\mum)';
% title('with volume filter')
xlabel('X (\mum)')
zlabel('Z (\mum)')
hold on
axis([-500 500 -125 125 -500 500])

figure(12)
quat_missing    = mapping_table_data.mapping_table(~idx_vol_diff,14:17)';
rod_missing     = ToFundamentalRegion(quat_missing, wscub.frmesh.symmetries);
% if ii == 1
PlotFRPerimeter('cubic');
% end
hold on
scatter3(rod_missing(1,:), rod_missing(2,:), rod_missing(3,:), ...
    20, mapping_table_data.mapping_table(~idx_vol_diff,22), 'filled')
axis equal tight off
title('Orientations of missing grains')
colormap jet
caxis([min(mapping_table_data.mapping_table(:,22))-0.1 max(mapping_table_data.mapping_table(:,22))+0.1])
hcb = colorbar;
set(hcb, 'Location', 'eastoutside')
hcb.Label.String    = 'grain radius (\mum)';

%%%
figure(2)
subplot(1,2,1)
plot(mapping_table_data.mapping_table(:,1), mapping_table_data.mapping_table(:,6), '.')
hold on
title('without volume filter')
xlabel('grain id (-)')
ylabel('misorientation (deg)')
axis tight

subplot(1,2,2)
plot(mapping_table_data.mapping_table(idx_vol_diff,1), mapping_table_data.mapping_table(idx_vol_diff,6), '.')
hold on
title('with volume filter')
xlabel('grain id (-)')
ylabel('misorientation (deg)')
axis tight

%%%
figure(3)
subplot(1,2,1)
plot(mapping_table_data.mapping_table(:,22), mapping_table_data.mapping_table(:,6), '.')
hold on
title('without volume filter')
xlabel('grain radius (\mum)')
ylabel('misorientation (deg)')
axis tight

subplot(1,2,2)
plot(mapping_table_data.mapping_table(idx_vol_diff,22), mapping_table_data.mapping_table(idx_vol_diff,6), '.')
hold on
title('with volume filter')
xlabel('grain radius (\mum)')
ylabel('misorientation (deg)')
axis tight

%%%
figure(4)
subplot(1,2,1)
plot(mapping_table_data.mapping_table(:,1), mapping_table_data.mapping_table(:,6)-median(mapping_table_data.mapping_table(:,6), 'omitnan'), '.')
hold on
title('without volume filter')
xlabel('grain id (-)')
ylabel('misorientation (deg)')
axis tight
% axis([0 400 -0.2 0.2])

subplot(1,2,2)
plot(mapping_table_data.mapping_table(idx_vol_diff,1), mapping_table_data.mapping_table(idx_vol_diff,6)-median(mapping_table_data.mapping_table(idx_vol_diff,6), 'omitnan'), '.')
hold on
title('with volume filter')
xlabel('grain id (-)')
ylabel('misorientation (deg)')
axis tight
% axis([0 400 -0.2 0.2])

%%%
figure(5)
subplot(1,2,1)
plot(mapping_table_data.mapping_table(:,22), mapping_table_data.mapping_table(:,23), '.')
title('without volume filter')
xlabel('grain radius from scan 0 (-)')
ylabel('grain radius from scan N (-)')
axis equal tight
hold on

subplot(1,2,2)
plot(mapping_table_data.mapping_table(idx_vol_diff,22), mapping_table_data.mapping_table(idx_vol_diff,23), '.')
title('with volume filter')
xlabel('grain radius from scan 0 (-)')
ylabel('grain radius from scan N (-)')
axis equal tight
hold on

figure, 
plot(mapping_table_data.mapping_table(:,3), mapping_table_data.mapping_table(:,6), 'o', ...
    'MarkerEdgeColor', 'r',...
    'MarkerFaceColor', 'r',...
    'MarkerSize', 5)
xlabel('grain number in S0 (-)')
ylabel('misorientation (deg)')
grid on

%%% DISTANCE
disp('dist')
disp([mean(mapping_table_data.mapping_table(:,5), 'omitnan') ...
    median(mapping_table_data.mapping_table(:,5), 'omitnan')  ...
    std(mapping_table_data.mapping_table(:,5), 'omitnan')])
disp([mean(mapping_table_data.mapping_table(idx_vol_diff,5), 'omitnan') ...
    median(mapping_table_data.mapping_table(idx_vol_diff,5), 'omitnan')  ...
    std(mapping_table_data.mapping_table(idx_vol_diff,5), 'omitnan')])

dist_stat               = [mean(mapping_table_data.mapping_table(:,5), 'omitnan') ...
    median(mapping_table_data.mapping_table(:,5), 'omitnan')  ...
    std(mapping_table_data.mapping_table(:,5), 'omitnan')];
dist_stat_vol_filter	= [mean(mapping_table_data.mapping_table(idx_vol_diff,5), 'omitnan') ...
    median(mapping_table_data.mapping_table(idx_vol_diff,5), 'omitnan')  ...
    std(mapping_table_data.mapping_table(idx_vol_diff,5), 'omitnan')];

%%% MISORIENTATION
disp('miso')
disp([mean(mapping_table_data.mapping_table(:,6), 'omitnan') ...
    median(mapping_table_data.mapping_table(:,6), 'omitnan')  ...
    std(mapping_table_data.mapping_table(:,6), 'omitnan')])
disp([mean(mapping_table_data.mapping_table(idx_vol_diff,6), 'omitnan') ...
    median(mapping_table_data.mapping_table(idx_vol_diff,6), 'omitnan')  ...
    std(mapping_table_data.mapping_table(idx_vol_diff,6), 'omitnan')])

miso_stat               = [mean(mapping_table_data.mapping_table(:,6), 'omitnan') ...
    median(mapping_table_data.mapping_table(:,6), 'omitnan')  ...
    std(mapping_table_data.mapping_table(:,6), 'omitnan')];
miso_stat_vol_filter	= [mean(mapping_table_data.mapping_table(idx_vol_diff,6), 'omitnan') ...
    median(mapping_table_data.mapping_table(idx_vol_diff,6), 'omitnan')  ...
    std(mapping_table_data.mapping_table(idx_vol_diff,6), 'omitnan')];

%%% X DISTANCE
disp('x-dist')
disp([mean(mapping_table_data.mapping_table(:,8) - mapping_table_data.mapping_table(:,11), 'omitnan') ...
    median(mapping_table_data.mapping_table(:,8) - mapping_table_data.mapping_table(:,11), 'omitnan') ...
    std(mapping_table_data.mapping_table(:,8) - mapping_table_data.mapping_table(:,11), 'omitnan')])
disp([mean(mapping_table_data.mapping_table(idx_vol_diff,8) - mapping_table_data.mapping_table(idx_vol_diff,11), 'omitnan') ...
    median(mapping_table_data.mapping_table(idx_vol_diff,8) - mapping_table_data.mapping_table(idx_vol_diff,11), 'omitnan') ...
    std(mapping_table_data.mapping_table(idx_vol_diff,8) - mapping_table_data.mapping_table(idx_vol_diff,11), 'omitnan')])

x_stat              = [mean(mapping_table_data.mapping_table(:,8) - mapping_table_data.mapping_table(:,11), 'omitnan') ...
    median(mapping_table_data.mapping_table(:,8) - mapping_table_data.mapping_table(:,11), 'omitnan') ...
    std(mapping_table_data.mapping_table(:,8) - mapping_table_data.mapping_table(:,11), 'omitnan')];
x_stat_vol_filter	= [mean(mapping_table_data.mapping_table(idx_vol_diff,8) - mapping_table_data.mapping_table(idx_vol_diff,11), 'omitnan') ...
    median(mapping_table_data.mapping_table(idx_vol_diff,8) - mapping_table_data.mapping_table(idx_vol_diff,11), 'omitnan') ...
    std(mapping_table_data.mapping_table(idx_vol_diff,8) - mapping_table_data.mapping_table(idx_vol_diff,11), 'omitnan')];

%%% Y DISTANCE
disp('y-dist')
disp([mean(mapping_table_data.mapping_table(:,9) - mapping_table_data.mapping_table(:,12), 'omitnan') ...
    median(mapping_table_data.mapping_table(:,9) - mapping_table_data.mapping_table(:,12), 'omitnan') ...
    std(mapping_table_data.mapping_table(:,9) - mapping_table_data.mapping_table(:,12), 'omitnan')])
disp([mean(mapping_table_data.mapping_table(idx_vol_diff,9) - mapping_table_data.mapping_table(idx_vol_diff,12), 'omitnan') ...
    median(mapping_table_data.mapping_table(idx_vol_diff,9) - mapping_table_data.mapping_table(idx_vol_diff,12), 'omitnan') ...
    std(mapping_table_data.mapping_table(idx_vol_diff,9) - mapping_table_data.mapping_table(idx_vol_diff,12), 'omitnan')])

y_stat              = [mean(mapping_table_data.mapping_table(:,9) - mapping_table_data.mapping_table(:,12), 'omitnan') ...
    median(mapping_table_data.mapping_table(:,9) - mapping_table_data.mapping_table(:,12), 'omitnan') ...
    std(mapping_table_data.mapping_table(:,9) - mapping_table_data.mapping_table(:,12), 'omitnan')];
y_stat_vol_filter	= [mean(mapping_table_data.mapping_table(idx_vol_diff,9) - mapping_table_data.mapping_table(idx_vol_diff,12), 'omitnan') ...
    median(mapping_table_data.mapping_table(idx_vol_diff,9) - mapping_table_data.mapping_table(idx_vol_diff,12), 'omitnan') ...
    std(mapping_table_data.mapping_table(idx_vol_diff,9) - mapping_table_data.mapping_table(idx_vol_diff,12), 'omitnan')];

disp('z-dist')
disp([mean(mapping_table_data.mapping_table(:,10) - mapping_table_data.mapping_table(:,13), 'omitnan') ...
    median(mapping_table_data.mapping_table(:,10) - mapping_table_data.mapping_table(:,13), 'omitnan') ...
    std(mapping_table_data.mapping_table(:,10) - mapping_table_data.mapping_table(:,13), 'omitnan')])
disp([mean(mapping_table_data.mapping_table(idx_vol_diff,10) - mapping_table_data.mapping_table(idx_vol_diff,13), 'omitnan') ...
    median(mapping_table_data.mapping_table(idx_vol_diff,10) - mapping_table_data.mapping_table(idx_vol_diff,13), 'omitnan') ...
    std(mapping_table_data.mapping_table(idx_vol_diff,10) - mapping_table_data.mapping_table(idx_vol_diff,13), 'omitnan')])

z_stat              = [mean(mapping_table_data.mapping_table(:,10) - mapping_table_data.mapping_table(:,13), 'omitnan') ...
    median(mapping_table_data.mapping_table(:,10) - mapping_table_data.mapping_table(:,13), 'omitnan') ...
    std(mapping_table_data.mapping_table(:,10) - mapping_table_data.mapping_table(:,13), 'omitnan')];
z_stat_vol_filter	= [mean(mapping_table_data.mapping_table(idx_vol_diff,10) - mapping_table_data.mapping_table(idx_vol_diff,13), 'omitnan') ...
    median(mapping_table_data.mapping_table(idx_vol_diff,10) - mapping_table_data.mapping_table(idx_vol_diff,13), 'omitnan') ...
    std(mapping_table_data.mapping_table(idx_vol_diff,10) - mapping_table_data.mapping_table(idx_vol_diff,13), 'omitnan')];

disp('median x-z-y disp')
disp([median(mapping_table_data.mapping_table(:,8) - mapping_table_data.mapping_table(:,11), 'omitnan') ...
    median(mapping_table_data.mapping_table(:,10) - mapping_table_data.mapping_table(:,13), 'omitnan') ...
    median(mapping_table_data.mapping_table(:,9) - mapping_table_data.mapping_table(:,12), 'omitnan')])
disp([median(mapping_table_data.mapping_table(idx_vol_diff,8) - mapping_table_data.mapping_table(idx_vol_diff,11), 'omitnan') ...
    median(mapping_table_data.mapping_table(idx_vol_diff,10) - mapping_table_data.mapping_table(idx_vol_diff,13), 'omitnan') ...
    median(mapping_table_data.mapping_table(idx_vol_diff,9) - mapping_table_data.mapping_table(idx_vol_diff,12), 'omitnan')])

xzy_median_stat             = [median(mapping_table_data.mapping_table(:,8) - mapping_table_data.mapping_table(:,11), 'omitnan') ...
    median(mapping_table_data.mapping_table(:,10) - mapping_table_data.mapping_table(:,13), 'omitnan') ...
    median(mapping_table_data.mapping_table(:,9) - mapping_table_data.mapping_table(:,12), 'omitnan')];
xzy_median_stat_vol_filter	= [median(mapping_table_data.mapping_table(idx_vol_diff,8) - mapping_table_data.mapping_table(idx_vol_diff,11), 'omitnan') ...
    median(mapping_table_data.mapping_table(idx_vol_diff,10) - mapping_table_data.mapping_table(idx_vol_diff,13), 'omitnan') ...
    median(mapping_table_data.mapping_table(idx_vol_diff,9) - mapping_table_data.mapping_table(idx_vol_diff,12), 'omitnan')];

disp('mean x-z-y disp')
disp([mean(mapping_table_data.mapping_table(:,8) - mapping_table_data.mapping_table(:,11), 'omitnan') ...
    mean(mapping_table_data.mapping_table(:,10) - mapping_table_data.mapping_table(:,13), 'omitnan') ...
    mean(mapping_table_data.mapping_table(:,9) - mapping_table_data.mapping_table(:,12), 'omitnan')])
disp([mean(mapping_table_data.mapping_table(idx_vol_diff,8) - mapping_table_data.mapping_table(idx_vol_diff,11), 'omitnan') ...
    mean(mapping_table_data.mapping_table(idx_vol_diff,10) - mapping_table_data.mapping_table(idx_vol_diff,13), 'omitnan') ...
    mean(mapping_table_data.mapping_table(idx_vol_diff,9) - mapping_table_data.mapping_table(idx_vol_diff,12), 'omitnan')])

xzy_mean_stat               = [mean(mapping_table_data.mapping_table(:,8) - mapping_table_data.mapping_table(:,11), 'omitnan') ...
    mean(mapping_table_data.mapping_table(:,10) - mapping_table_data.mapping_table(:,13), 'omitnan') ...
    mean(mapping_table_data.mapping_table(:,9) - mapping_table_data.mapping_table(:,12), 'omitnan')];
xzy_mean_stat_vol_filter	= [mean(mapping_table_data.mapping_table(idx_vol_diff,8) - mapping_table_data.mapping_table(idx_vol_diff,11), 'omitnan') ...
    mean(mapping_table_data.mapping_table(idx_vol_diff,10) - mapping_table_data.mapping_table(idx_vol_diff,13), 'omitnan') ...
    mean(mapping_table_data.mapping_table(idx_vol_diff,9) - mapping_table_data.mapping_table(idx_vol_diff,12), 'omitnan')];