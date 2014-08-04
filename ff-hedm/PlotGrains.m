clear all
close all
clc

pfname  = 'O:\Sharma_Oct13\Stebner_Dec13\Ring4\Grains.csv';
grains  = load(pfname);

idx = (abs(grains(:, 13)) < 180);

figure(1) 
scatter3(grains(:,11), grains(:,12), grains(:,13), 10, 'filled') %% COMPLETENESS
figure(2) 
scatter3(grains(idx,11), grains(idx,12), grains(idx,13), 10, 'filled')
return
pname0  = 'O:\stebner_dec13\wedge.tx.map';

wedge   = -1.0:0.5:1.0;
tx      = -1.0:0.5:1.0;

for i = 3%1:1:length(wedge)
    for j = 1:1:length(tx)
        pname1  = ['wedge_', num2str(wedge(i), '%1.1f'), '_tx_', num2str(tx(i), '%1.1f')];
        pname2  = 'PeakSearch\NDC5_1';
        fname   = 'Grains.csv';
        pfname  = fullfile(pname0, pname1, pname2, fname);
        grains  = load(pfname);
        idx = (abs(grains(:, 13)) < 180);
        
        figure(i*10+j)
        % scatter3(grains(:,11), grains(:,12), grains(:,13), 10, grains(:,24), 'filled') %% COMPLETENESS
        scatter3(grains(idx,11), grains(idx,12), grains(idx,13), 10, grains(idx,24), 'filled') %% COMPLETENESS
        pause
    end
end

return
%%%%%%% : ring1 : ring2 : ring3
% ring1 : 501   :   285 : 
% ring2 :       :   596 :
% ring3 :       :       :   181 

qsym    = CubSymmetries;

pfname  = '.\Grains_rings1_2_3_4_5.csv';
grains0 = load(pfname);

nGrains0    = size(grains0, 1);

for i = 1:1:nGrains0
    RMats0(:,:,i)   =  reshape(grains0(i,2:10), 3, 3)';
end
quat0   = ToFundamentalRegionQ(QuatOfRMat(RMats0), qsym);
rod0    = RodOfQuat(quat0);
COM0    = grains0(:,11:13);

%%%% PLOTTING
figure(1)
% scatter3(grains0(:,11), grains0(:,12), grains0(:,13), 10,
% grains0(:,20)./max(grains0(:,20)), 'filled') %% ERROR
scatter3(grains0(:,11), grains0(:,12), grains0(:,13), 10, grains0(:,24), 'filled') %% COMPLETENESS

figure(2)
PlotFRPerimeter('cubic');
plot3(rod0(1,:), rod0(2,:), rod0(3,:), 'k.')
return

% Orientation: O[row][col]
% Sp_ID O[0][0] O[0][1] O[0][2] O[1][0] O[1][1] O[1][2] O[2][0] O[2][1] O[2][2] X Y Z a b c alpha beta gamma Err1 Err2 Err3 MeanRadius Confidence

qsym    = CubSymmetries;

pname   = 'O:\obstalecki_may14\parkjs\Ring1';
fname0  = 'PeakSearch\sample_1_7_1\Grains.csv';
pfname0 = fullfile(pname, fname0);

pname   = 'O:\obstalecki_may14\parkjs\Ring3';
fname1  = 'PeakSearch\sample_1_7_1\Grains.csv';
pfname1 = fullfile(pname, fname1);

grains0 = load(pfname0);
grains1 = load(pfname1);

% grains0 = grains0(1:12, :);
% grains1 = grains1(1:8, :);

nGrains0    = size(grains0, 1);
nGrains1    = size(grains1, 1);

for i = 1:1:nGrains0
    RMats0(:,:,i)   =  reshape(grains0(i,2:10), 3, 3)';
end
quat0   = ToFundamentalRegionQ(QuatOfRMat(RMats0), qsym);
rod0    = RodOfQuat(quat0);
COM0    = grains0(:,11:13);

for i = 1:1:nGrains1
    RMats1(:,:,i)   =  reshape(grains1(i,2:10), 3, 3)';
end
quat1   = ToFundamentalRegionQ(QuatOfRMat(RMats1), qsym);
rod1    = RodOfQuat(quat1);
COM1    = grains1(:,11:13);

MisoAngle       = zeros(nGrains0, nGrains1);
MisoAngleThresh = 0.5;
for i = 1:1:nGrains1
    MisoAngle(:,i)  = Misorientation(quat1(:,i), quat0, qsym);
end
MisoAngle       = rad2deg(MisoAngle);
[MisoAngleMin, MisoAngleMinIdx] = min(MisoAngle, [], 2);

dCOM        = zeros(nGrains0, nGrains1);
dCOMThresh  = 100;
for i = 1:1:nGrains1
    dCOMi   = (COM0 - repmat(COM1(i,:), nGrains0, 1)).^2;
    dCOM(:,i)   = sqrt(sum(dCOMi,2));
end
[dCOMMin, dCOMMinIdx]   = min(dCOM, [], 2);

Flag_MisoAngle  = find(MisoAngle < MisoAngleThresh);
Flag_dCOM       = find(dCOM < dCOMThresh);
Flag            = (MisoAngle < MisoAngleThresh) & (dCOM < dCOMThresh);

%%%% PLOTTING
figure(1)
subplot(1,2,1)
scatter3(grains0(:,11), grains0(:,12), grains0(:,13), 10, grains0(:,20)./max(grains0(:,20)), 'filled')
subplot(1,2,2)
scatter3(grains1(:,11), grains1(:,12), grains1(:,13), 10, grains1(:,20)./max(grains1(:,20)), 'filled')

figure(2)
subplot(1,2,1)
PlotFRPerimeter('cubic');
axis square off
plot3(rod0(1,:), rod0(2,:), rod0(3,:), 'k.')

subplot(1,2,2)
PlotFRPerimeter('cubic');
axis square off
plot3(rod1(1,:), rod1(2,:), rod1(3,:), 'k.')

% figure(3)
% subplot(1,2,1)
% scatter3(grains0(:,11), grains0(:,12), grains0(:,13), 10, grains0(:,20)./max(grains0(:,20)), 'filled')