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

%%% GENERATE SH MODES
if InputData.UseProjectedHarmonics
    disp('projecting standard harmonics to frmesh ...')
    load wshex12x
    dh  = zeros(frmesh.numind, InputData.nHarm);
    for i = 1:1:InputData.nHarm
        dh(:,i) = Map2StdMesh( ...
            wshex12x.frmesh.dh(:,i), ...
            wshex12x.frmesh, ...
            frmesh);
    end
    
    clear wshex12x wsopts
    disp('done')
else
    [dh, ~]     = DiscreteHarmonics(frmesh, InputData.nHarm);
end

%%% ASSEMBLE DH MATRIX
DH	= zeros(6*frmesh.numind,6*InputData.nHarm);

DH(1:6:end-5,1:InputData.nHarm)                     = dh;           % 11
DH(2:6:end-4,InputData.nHarm+1:2*InputData.nHarm)   = dh;           % 22
DH(3:6:end-3,2*InputData.nHarm+1:3*InputData.nHarm) = dh;           % 33
DH(4:6:end-2,3*InputData.nHarm+1:4*InputData.nHarm) = dh;           % 12
DH(5:6:end-1,4*InputData.nHarm+1:5*InputData.nHarm) = dh;           % 13
DH(6:6:end,5*InputData.nHarm+1:6*InputData.nHarm)   = dh;           % 23

DH	= sparse(DH);

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
    Sigma_d = DH*w_di;
    Sigma_d = Sigma_d;     
            
    ig          = IG(ri:rf).*1000000;           %%% TPa <> MPa
    Sigma_d_ig  = DH*ig;
    Sigma_d_ig  = Sigma_d_ig;
    
    ub          = UB(ri:rf).*1000000;
    lb          = LB(ri:rf).*1000000;
    
    Sigma_d     = reshape(Sigma_d, 6, frmesh.numind);
    Sigma_d_ig  = reshape(Sigma_d_ig, 6, frmesh.numind);
    
    w_di    = reshape(w_di, InputData.nHarm, 6);
    ig      = reshape(ig, InputData.nHarm, 6);
    ub      = reshape(ub, InputData.nHarm, 6);
    lb      = reshape(lb, InputData.nHarm, 6);
    
    Sigma_d(4:6,:)      = Sigma_d(4:6,:)./sqrt(2);
    Sigma_d_ig(4:6,:)   = Sigma_d_ig(4:6,:)./sqrt(2);
    for j = 1:1:6
        PlotFR(frmesh, Sigma_d(j,:), ...
            'Symmetries', 'hexagonal', ...
            'ShowMesh', 'on')
        PlotFR(frmesh, Sigma_d_ig(j,:), ...
            'Symmetries', 'hexagonal', ...
            'ShowMesh', 'on')
    end
end