clear all
close all
clc

wsname  = 'wscub4x';      % workspace name

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECKING REPEATS
% % %%% LAYER A
% pname1  = '/home/beams/S1IDUSER/mnt/orthros/gao_aug15/hedm_correct_a/RP1A_Layer1_Analysis_Time_2016_10_13_13_09_16';
% fname1  = 'Grains.csv';
% % %%% LAYER B
% pname2  = '/home/beams/S1IDUSER/mnt/orthros/gao_aug15/hedm_correct_a/RP1A_Layer5_Analysis_Time_2016_10_14_11_31_51';
% fname2  = 'Grains.csv';
%%%
% %%% LAYER A
pname1  = '/home/beams/S1IDUSER/mnt/orthros/gao_aug15/hedm_correct_a/RP4A_Layer1_Analysis_Time_2016_10_14_21_14_50';
fname1  = 'Grains.csv';
% %%% LAYER B
pname2  = '/home/beams/S1IDUSER/mnt/orthros/gao_aug15/hedm_correct_a/RP4A_Layer5_Analysis_Time_2016_10_16_21_21_36';
fname2  = 'Grains.csv';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RLab2Sam    = eye(3,3);

% COLORING SCHEME
% ColorMap = load('cubic_xstal_coloring_scheme.mat');

% XSTAL SYMMETRY IN MTEX CONVENTION
% cs  = symmetry('m-3m');

% SAMPLE SYMMETRY IN MTEX CONVENTION
% ss  = symmetry('-1');

% FILTERS
Thresh_Completeness = 0.8;
Thresh_GrainRadius  = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Load workspace for fundamental region.
load(wsname);
eval(['ws = ', wsname, ';']);
clear(wsname)

% Load MIDAS results
pfname1 = fullfile(pname1, fname1);
Grains1 = parseGrainData(pfname1, ws.frmesh.symmetries);

% Load MIDAS results
pfname2 = fullfile(pname2, fname2);
Grains2 = parseGrainData(pfname2, ws.frmesh.symmetries);

% THRESHOLDING BY COMPLETENESS
idx_Completeness1   = [Grains1.Completeness] >= Thresh_Completeness;
idx_MeanRadius1     = [Grains1.GrainRadius] >= Thresh_GrainRadius;
idx1    = find(idx_Completeness1);
numpts1 = length(idx1);
wts1    = ones(1, numpts1);

grainID1    = [Grains1(idx1).GrainID]';
xyz1        = [Grains1(idx1).COM]';
rod1        = [Grains1(idx1).rod];
cidx1       = [Grains1(idx1).Completeness];
quat1       = [Grains1(idx1).quat];
quat1       = ToFundamentalRegionQ(quat1, ws.frmesh.symmetries);
GrainRad1   = [Grains1(idx1).GrainRadius];
lattprm1    = [Grains1(idx1).lattprms];
compl1      = [Grains1(idx1).Completeness];

% THRESHOLDING BY COMPLETENESS
idx_Completeness2   = [Grains2.Completeness] >= Thresh_Completeness;
idx_MeanRadius2     = [Grains2.GrainRadius] >= Thresh_GrainRadius;
idx2    = find(idx_Completeness2);
numpts2 = length(idx2);
wts2    = ones(1, numpts2);

grainID2    = [Grains2(idx2).GrainID]';
xyz2        = [Grains2(idx2).COM]';
rod2        = [Grains2(idx2).rod];
cidx2       = [Grains2(idx2).Completeness];
quat2       = [Grains2(idx2).quat];
quat2       = ToFundamentalRegionQ(quat2, ws.frmesh.symmetries);
GrainRad2   = [Grains2(idx2).GrainRadius];
lattprm2    = [Grains2(idx2).lattprms];
compl2      = [Grains2(idx2).Completeness];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOOK FOR GRAINS WITH SIMILAR ORIENTATION & COM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xyz2_temp   = xyz2;
quat2_temp  = quat2;
DiffChart   = cell(numpts1,1);
DiffTable   = zeros(numpts1,8);
for i = 1:1:numpts1
    xyz1i   = xyz1(i,:);
    rod1i   = rod1(:,i);
    quat1i  = quat1(:,i);
    
    distance    = sqrt((xyz2(:,1) - xyz1i(1)).^2 + (xyz2(:,2) - xyz1i(2)).^2 + (xyz2(:,3) - xyz1i(3)).^2);
    miso        = Misorientation(quat1i, quat2, ws.frmesh.symmetries)';
    miso        = rad2deg(miso);
    
    [distance_sorted, idx_distance] = sort(distance, 'ascend');
    [miso_sorted, idx_miso]         = sort(miso, 'ascend');
    
    DiffChart{i,1}  = [idx_distance, distance_sorted, idx_miso, miso_sorted];
    
    [dist_min, dist_idx]    = min(distance);
    [miso_min, miso_idx]   	= min(miso);
    dist_min_based_on_miso  = miso(miso_idx);
    
%     xyz2_temp(dist_idx, :)  = nan;
%     quat2_temp(:, mo_idx)   = nan;
    
    DiffTable(i,1)  = i;
    DiffTable(i,2)  = dist_min;
    DiffTable(i,3)  = dist_idx;
    DiffTable(i,4)  = miso_min;
    DiffTable(i,5)  = miso_idx;
    DiffTable(i,6)  = dist_min_based_on_miso;
    DiffTable(i,7)  = dist_min_based_on_miso < 200;
    DiffTable(i,8)  = compl1(i);
end

idx_similar = find(DiffTable(:,3) == DiffTable(:,5));
idx_diff    = find(DiffTable(:,3) ~= DiffTable(:,5));

% for i = 1:1:length(idx_diff)
%     % DiffChart{idx_diff(i),1}(1:10,:)
%     [DiffChart{idx_diff(i),1}(1:10,1) DiffChart{idx_diff(i),1}(1:10,3)]
% 	pause
% end
figure,
plot(GrainRad1(idx_similar), 'b.')
hold on
plot(GrainRad1(idx_diff), 'r.')
xlabel('grain counter')
ylabel('grain radius (um)')
legend('similar', 'dissimilar')

figure, 
PlotFRPerimeter('hexagonal')
axis equal
hold on
scatter3(rod1(1,idx_similar), rod1(2,idx_similar), rod1(3,idx_similar), 10, 'b', 'filled')
scatter3(rod1(1,idx_diff), rod1(2,idx_diff), rod1(3,idx_diff), 10, 'r', 'filled')
% scatter3(rod1(1,idx_similar), rod1(2,idx_similar), rod1(3,idx_similar), 30, DiffTable(idx_similar,6), 'filled')
% scatter3(rod1(1,idx_diff), rod1(2,idx_diff), rod1(3,idx_diff), 30, DiffTable(idx_diff,6), 'filled')

figure, 
scatter3(xyz1(idx_similar,1), xyz1(idx_similar,2), xyz1(idx_similar,3), 10, 'b', 'filled')
hold on
scatter3(xyz1(idx_diff,1), xyz1(idx_diff,2), xyz1(idx_diff,3), 10, 'r', 'filled')
% scatter3(xyz1(idx_similar,1), xyz1(idx_similar,2), xyz1(idx_similar,3), 30, DiffTable(idx_similar,6), 'filled')
% hold on
% scatter3(xyz1(idx_diff,1), xyz1(idx_diff,2), xyz1(idx_diff,3), 30, DiffTable(idx_diff,6), 'filled')
axis square
legend('similar', 'dissimilar')

disp(sprintf('grains in the first set  : %d', numpts1));
disp(sprintf('grains in the second set : %d', numpts2));
disp(sprintf('A. similar grains based on miso and COM              : %d', length(idx_similar)));
disp(sprintf('B. different grains based on miso and COM            : %d', length(idx_diff)));
disp(sprintf('C. different grains based on miso and (COM < 0.2 mm) : %d', sum(DiffTable(:,7))));
disp(sprintf('mean / std dist diff between similar grains (A) : %f %f', [mean(DiffTable(idx_similar,2)), std(DiffTable(idx_similar,2))]))
disp(sprintf('mean / std ang diff between similar grains  (A) : %f %f', rad2deg([mean(DiffTable(idx_similar,4)), std(DiffTable(idx_similar,4))])))
disp(sprintf('mean / std grain radius of similar grains   (A) : %f %f', [mean(GrainRad1(idx_similar)), std(GrainRad1(idx_similar))]))
disp(sprintf('mean / std dist diff between similar grains (C) : %f %f', [mean(DiffTable(DiffTable(:,7),2)), std(DiffTable(DiffTable(:,7),2))]))
disp(sprintf('mean / std ang diff between similar grains  (C) : %f %f', rad2deg([mean(DiffTable(DiffTable(:,7),4)), std(DiffTable(DiffTable(:,7),4))])))
disp(sprintf('mean / std grain radius of similar grains   (C) : %f %f', [mean(GrainRad1(DiffTable(:,7))), std(GrainRad1(DiffTable(:,7)))]))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HISTOGRAM OF GRAIN SIZES NORMALIZED BY MAX GRAIN SIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbins   = 20;
binsab    = linspace(0, 1, nbins);
h11      = hist(GrainRad1./max(GrainRad1), binsab)';
h21      = hist(GrainRad2./max(GrainRad2), binsab)';

p0  = [150; 0.2; 0.2; 0.5; 0.5; 0.3];
h1fit0  = pksplitpseudoVoigt(p0, binsab);
p1      = lsqcurvefit(@pksplitpseudoVoigt, p0, binsab, h11);
h1fit   = pksplitpseudoVoigt(p1, binsab);

h2fit0  = pksplitpseudoVoigt(p0, binsab);
p2      = lsqcurvefit(@pksplitpseudoVoigt, p0, binsab, h21);
h2fit   = pksplitpseudoVoigt(p2, binsab);

figure,
bar(binsab, [h11 h21]);
hold on
% plot(bins, h1fit0, 'b:')
plot(binsab, h1fit, 'b-')
% plot(bins, h2fit0, 'r:')
plot(binsab, h2fit, 'r-')
xlabel('relative grain radius (-)')
ylabel('number of grains (-)')
grid on
title(sprintf('Max grain size : coarse scan %3.0f (um) | fine scan %3.0f (um)', max(GrainRad1), max(GrainRad2)))
legend('load 4', 'load N', 'coarse scan fit', 'fine scan fit')
axis([0 1 0 300])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HISTOGRAM OF LATTICE PARAMETER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbins   = 20;
ab1     = [lattprm1(1,:) lattprm1(2,:)]';
c1      = lattprm1(3,:)';
ab2     = [lattprm2(1,:) lattprm2(2,:)]';
c2      = lattprm2(3,:)';
binsab  = linspace(min([ab1; ab2]), max([ab1; ab2]), nbins);
binsc   = linspace(min([c1; c2]), max([c1; c2]), nbins);
hab1    = hist(ab1, binsab)';
hab2    = hist(ab2, binsab)';
hc1     = hist(c1, binsc)';
hc2     = hist(c2, binsc)';

figure,
bar(binsab, [hab1 hab2])
grid on
xlabel('a / b (Angstrom)')
ylabel('number of grains (-)')
legend('load 4', 'load N')

figure,
bar(binsab, [hc1 hc2])
grid on
xlabel('c (Angstrom)')
ylabel('number of grains (-)')
legend('load 4', 'load N')

nbins   = 20;
angab1  = [lattprm1(4,:) lattprm1(5,:)]';
angc1   = lattprm1(6,:)';
angab2  = [lattprm2(4,:) lattprm2(5,:)]';
angc2   = lattprm2(6,:)';
binsab  = linspace(min([angab1; angab2]), max([angab1; angab2]), nbins);
binsc   = linspace(min([angc1; angc2]), max([angc1; angc2]), nbins);
hab1    = hist(angab1, binsab)';
hab2    = hist(angab2, binsab)';
hc1     = hist(angc1, binsc)';
hc2     = hist(angc2, binsc)';

figure,
bar(binsab, [hab1 hab2])
grid on
xlabel('alpha / beta (deg)')
ylabel('number of grains (-)')
legend('load 4', 'load N')

figure,
bar(binsab, [hc1 hc2])
grid on
xlabel('gamma (deg)')
ylabel('number of grains (-)')
legend('load 4', 'load N')