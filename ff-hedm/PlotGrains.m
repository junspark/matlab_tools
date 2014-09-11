clear all
close all
clc

%%% FIRST RUN : startup_mtex.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%v
%%% GENERATE IPF COLORMAP USING MTEX
%%% CHECK IN ebsdColorbar.m
% cs  = symmetry('m-3m');
% ss  = symmetry('-1');
% cc  = get_option('antipodal','colorcoding','ipdfHKL');

% [minTheta,maxTheta,minRho,maxRho,v] = getFundamentalRegionPF(cs, 'antipodal');
% h   = S2Grid('PLOT', 'minTheta', minTheta, 'maxTheta', maxTheta,...
%     'minRho', minRho, 'maxRho', maxRho, 'RESTRICT2MINMAX', 'resolution', 1*degree, 'antipodal');
% v   = vector3d(h);
% x   = getx(v); x = x(:);
% y   = gety(v); y = y(:);
% z   = getz(v); z = z(:);
% d   = orientation2color(h,cc,cs,'antipodal');
% save('coloring_scheme.mat', 'x', 'y', 'z', 'd')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%v

load('.\coloring_scheme.mat');
% create an EBSD variable containing the data
% DATA
% Sp_ID O[0][0] O[0][1] O[0][2] O[1][0] O[1][1] O[1][2] O[2][0] O[2][1] O[2][2] X Y Z a b c alpha beta gamma Err1 Err2 Err3 MeanRadius Confidence 
grains	= load('.\grains_example.csv');

nGrains     = size(grains, 1);
for i = 1:1:nGrains
    RMats(:,:,i)   =  reshape(grains(i,2:10), 3, 3)';
end
qsym    = CubSymmetries; Rsym    = RMatOfQuat(qsym);
quat    = ToFundamentalRegionQ(QuatOfRMat(RMats), qsym);
rod     = RodOfQuat(quat);

ori     = orientation('quaternion', quat(1,:), quat(2,:), quat(3,:), quat(3,:), cs);
ebsd    = EBSD(ori, cs, ss);
rgb     = orientation2color(ori, 'ipdfHSV');

plot(ebsd, 'colorcoding','ipdfHSV')

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

%%%% PLOT COM / RGB IN FUNDAMENTAL TRIANGLE AS IPDF
figure, scatter3(grains(:,11), grains(:,12), grains(:,13), 50, rgb, 'filled') %% COMPLETENESS
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

%%%% PLOT ORIENTATIONS / RGB IN FUNDAMENTAL TRIANGLE AS IPDF
figure, PlotFRPerimeter('cubic');
scatter3(rod(1,:), rod(2,:), rod(3,:), 50, rgb, 'filled')
axis square tight off