clear all
close all
clc

ndiv    = 1000;

logname     = '/home/jspark/Desktop/fablework/OMC_S1_00013/peaksearch_OMC_S1_00013.log';
GrainData   = parseGrainSpotterLog(logname);
num_grain	= length(GrainData);

%%% SPECIAL DIRECTIONS
LD      = [0 1 0]';
stiff_fcc   = [1 1 1]';
compl_fcc   = [1 0 0]';     %%% NEED TO CORRECT

%%% GENERATE FIBERS ONCE & LOAD 
qsym_fcc    = CubSymmetries;
fib_stiff	= QuatOfRod(FiberOfPoint(LD, stiff_fcc, ndiv, qsym_fcc));
fib_stiff	= ToFundamentalRegionQ(fib_stiff, qsym_fcc);
fib_stiff_rod   = RodOfQuat(fib_stiff);

fib_compl	= QuatOfRod(FiberOfPoint(LD, compl_fcc, ndiv, qsym_fcc));
fib_compl	= ToFundamentalRegionQ(fib_compl, qsym_fcc);
fib_compl_rod   = RodOfQuat(fib_compl);

subplot(1,2,1)
PlotFRPerimeter('cubic');
hold on
axis equal tight
plot3(fib_stiff_rod(1,:), fib_stiff_rod(2,:), fib_stiff_rod(3,:), 'k.')

subplot(1,2,2)
PlotFRPerimeter('cubic');
hold on
axis equal tight
plot3(fib_compl_rod(1,:), fib_compl_rod(2,:), fib_compl_rod(3,:), 'k.')

angle_stiff = zeros(num_grain,1);
angle_compl = zeros(num_grain,1);
for i = 1:1:num_grain
    GrainData(i).completeness  = GrainData(i).nMeasGvec./GrainData(i).nExpGvec; 
    
    quats_gr	= ToFundamentalRegionQ(QuatOfRMat(GrainData(i).U), qsym_fcc);    %%% CHECK TRANSPOSE OR NOT U IN FABLE CONVENTION
    
    anglei  = Misorientation(quats_gr, fib_stiff, qsym_fcc);
    anglei  = min(anglei);
    angle_stiff(i,1)	= anglei;
    
    anglei  = Misorientation(quats_gr, fib_compl, qsym_fcc);
    anglei  = min(anglei);
    angle_compl(i)  = anglei;
end
GrainData(i).angle_stiff	= angle_stiff;      %%% MISORIENTATION WRT STIFF || LD FIB
GrainData(i).angle_compl	= angle_compl;      %%% MISORIENTATION WRT COMPLIANT || LD FIB

idx_stiff           = find(rad2deg(GrainData(i).angle_stiff) < 10);
stiff               = angle_stiff(idx_stiff);
stiff_n_complete    = [GrainData(idx_stiff).completeness]';
[idx_stiff rad2deg(stiff) stiff_n_complete]

idx_compl           = find(rad2deg(GrainData(i).angle_compl) < 10);
compl               = angle_compl(idx_compl);
compl_n_complete    = [GrainData(idx_compl).completeness]';
[idx_compl rad2deg(compl) compl_n_complete]