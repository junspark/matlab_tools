clear all
close all
clc

% Sp_ID O[0][0] O[0][1] O[0][2] O[1][0] O[1][1] O[1][2] O[2][0] O[2][1] O[2][2] X Y Z a b c alpha beta gamma Err1 Err2 Err3 MeanRadius Confidence 
qsym    = CubSymmetries;
grains	= load('.\grains_example.csv');

nGrains    = size(grains, 1);
for i = 1:1:nGrains
    RMats(:,:,i)   =  reshape(grains(i,2:10), 3, 3)';
end
quat    = ToFundamentalRegionQ(QuatOfRMat(RMats), qsym);
rod     = RodOfQuat(quat);


%%% THRESHOLDING BY COMPLETENESS
Thresh_Completeness = 0.7;
idx_Completeness    = grains(:,24) >= Thresh_Completeness;

%%% THRESHOLDING BY MEAN RADIUS
Thresh_MeanRadius   = 50;
idx_MeanRadius      = grains(:,23) >= Thresh_MeanRadius;

%%%% PLOT COM / ONE COLOR
figure, scatter3(grains(:,11), grains(:,12), grains(:,13), 50, 'filled', 'b')
grid on
axis square

%%%% PLOT COM / COMPLETENESS AS COLOR
figure, scatter3(grains(:,11), grains(:,12), grains(:,13), 50, grains(:,24), 'filled') %% COMPLETENESS
grid on
axis square

%%%% PLOT ORIENTATIONS / ONE COLOR
figure, PlotFRPerimeter('cubic');
scatter3(rod(1,:), rod(2,:), rod(3,:), 50, 'filled', 'b')
axis square tight off

%%%% PLOT ORIENTATIONS / COMPLETENESS AS COLOR
figure, PlotFRPerimeter('cubic');
scatter3(rod(1,:), rod(2,:), rod(3,:), 50, grains(:,24), 'filled')
axis square tight off
colorbar vert