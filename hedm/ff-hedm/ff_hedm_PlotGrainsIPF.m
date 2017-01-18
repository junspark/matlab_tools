clear all
close all
clc

%%% INSTALLS MTEX
run('/home/beams/S1IDUSER/mnt/s1b/__eval/mtex-4.4.0/startup_mtex.m')
data    = load('grains_unirr_def.csv');

cs  = crystalSymmetry('Td');

R_ESRF2APS  = RMatOfQuat(QuatOfESRF2APS);

cluster_id  = data(:,1);
cid         = 210;

idx     = find(cluster_id == cid);
RMATs   = data(idx, 2:10);

figure(10) 
PlotFRPerimeter('cubic')
axis square off
hold on
for j = 1:1:length(idx)
    rmat_c2l_ESRF	= reshape(RMATs(j,:), 3, 3)';
    rmat_c2l_APS    = R_ESRF2APS*rmat_c2l_ESRF;
    rmat_l2c_APS    = rmat_c2l_APS';
    
    % CHECK SAMPLE TO XSTAL CONVENTION
    % IN CONVENTIONAL TEXTURE ANALYSIS ROTATION IS DEFINED AS Rs = c
    Quat    = ToFundamentalRegionQ(QuatOfRMat(rmat_l2c_APS), CubSymmetries);
    R       = RMatOfQuat(Quat);
    bunge(:,j)  = BungeOfRMat(R, 'radians');
    
    %%% PLOT ORIENTATION IN RODRIGUES SPACE
    Rod     = RodOfQuat(ToFundamentalRegionQ(QuatOfRMat(rmat_c2l_APS), CubSymmetries));
    figure(10)
    scatter3(Rod(1), Rod(2), Rod(3), 'k.')
end
figure(10)
hold off

ori = orientation('Euler', bunge(1,:), bunge(2,:), bunge(3,:), cs);

figure,
plotIPDF(ori, yvector, cs, 'antipodal')
hold off

%%% UNINSTALLS MTEX IN CASE OF CONFLICTS
run('/home/beams/S1IDUSER/mnt/s1b/__eval/mtex-4.4.0/uninstall_mtex.m')