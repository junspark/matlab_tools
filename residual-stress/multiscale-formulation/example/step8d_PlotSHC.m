clear all
close all
clc

%%% INPUT DATA GENERATED IN STEP0
InputData   = load('example_rsinput.mat');
frmesh      = InputData.frmesh;

%%% FILE NAME FOR SOLUTION
FNAME_SOL   = 'Solution_l1_1_l2_0.1_l3_0.1_k_1.5.mat';
PFNAME_SOL  = fullfile(InputData.PNAME_SOLUTION, FNAME_SOL);
SOL = load(PFNAME_SOL);
SOL = SOL.SOL2;

%%% INITIAL GUESS
PFNAME_INVIG    = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_INVIG);
load(PFNAME_INVIG)

%%% DV INFO
PFNAME_GRID = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_GRID);
GridData    = load(PFNAME_GRID);

%%% MESH INFO
PFNAME_MESH = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_MESH);
MeshData    = load(PFNAME_MESH);

Sigma_h = SOL(1:6*MeshData.numnp);
w_d     = SOL(6*MeshData.numnp+1:end);
for i = 1:1:GridData.ndv
    odf_data    = load(InputData.PFNAME_ODF{i});
    PlotFR(frmesh, odf_data, ...
        'Symmetries', 'hexagonal', ...
        'ShowMesh', 'on')
    
    ri  = 1 + 6*(i-1)*InputData.nHarm;
    rf  = 6*i*InputData.nHarm;
    
    w_di    = w_d(ri:rf).*1000000;              %%% TPa <> MPa
    ig      = IG(ri:rf).*1000000;           %%% TPa <> MPa
    ub      = UB(ri:rf).*1000000;           %%% TPa <> MPa
    lb      = LB(ri:rf).*1000000;           %%% TPa <> MPa
    
    %%% THESE ARE STILL IN XX,YY,ZZ,YZ,XZ,XY ORDER
    w_di    = reshape(w_di, InputData.nHarm, 6);
    ig      = reshape(ig, InputData.nHarm, 6);
    ub      = reshape(ub, InputData.nHarm, 6);
    lb      = reshape(lb, InputData.nHarm, 6);
    
    w_di(:,4:6) = w_di(:,4:6)./sqrt(2);
    ig(:,4:6)   = ig(:,4:6)./sqrt(2);
    ub(:,4:6)   = ub(:,4:6)./sqrt(2);
    lb(:,4:6)	= lb(:,4:6)./sqrt(2);
    
    for j = 1:1:6
        figure,
        plot(1:InputData.nHarm, ig(:,j), 'rs-')
        hold on
        plot(1:InputData.nHarm, ub(:,j), 'm^-')
        plot(1:InputData.nHarm, lb(:,j), 'mv-')
        plot(1:InputData.nHarm, w_di(:,j), 'go-')
        axis tight
        xlabel('SH mode number')
        ylabel('SH coefficient (MPa)')
        legend({'IG'; 'UB'; 'LB'; 'w'}, 'location', 'best')
    end
    return
    close all
end