clear all
close all
clc

%%% SOLVER COMBINING K Q A AND N
% RUN ONCE AND SAVE RESULT OR ELSE WASTE HUGE TIME
% MIGHT HAVE TO MIGRATE TO BIGGER SOLVER

%%% NEED TO APPLY TRACTIONS FOR TOP AND BOTTOM FACES

%%% INPUT DATA GENERATED IN STEP0
InputData   = load('example_rsinput.mat');

frmesh  = InputData.frmesh;

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

%%% FILE NAME FOR SOLUTION
FNAME_SOL   = 'Solution_l1_1_l2_0.1_l3_0.1_k_1.5.mat';
PFNAME_SOL  = fullfile(InputData.PNAME_SOLUTION, FNAME_SOL);

%%% MESH INFO
PFNAME_MESH = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_MESH);
MeshData    = load(PFNAME_MESH);

%%% DV INFO
PFNAME_GRID = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_GRID);
GridData    = load(PFNAME_GRID);

%%% MICRO MATRIX
PFNAME_A    = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_A);
load(PFNAME_A)
PFNAME_N    = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_N);
load(PFNAME_N)

%%% MACRO MATRICES
PFNAME_KQ    = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_KQ);
load(PFNAME_KQ)

PFNAME_Ksym = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_Ksym);
load(PFNAME_Ksym)

PFNAME_Kfs  = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_Kfs);
load(PFNAME_Kfs)

PFNAME_Keq  = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_Keq);
load(PFNAME_Keq)

%%% VECTORS
PFNAME_Fvec     = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_Fvec);
load(PFNAME_Fvec)
PFNAME_INVIG    = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_INVIG);
load(PFNAME_INVIG)

%%% OBTAIN SIZES
[nA, mA]    = size(A);
[nK, mK]    = size(K);
[nN, mN]    = size(N);
[nF, mF]    = size(Fvec);
[nIG, mIG]  = size(IG);

%%% METHOD 1 - SOLVE SIMULTANEOUSLY
Fmacro1 = repmat(InputData.SIMGA_APPLIED, 1, MeshData.numnp);
Fmacro1 = Fmacro1(:);

F       = [Fmacro1; Fvec];

LHS     = [ ...
    K             -Q*N; ...
    zeros(nA, mK)    A; ...
    ];

SOL1    = LHS\F;

SIGMA1_h    = SOL1(1:nK,1);
SH1         = SOL1(nK+1:end,:);

disp(norm(Keq*SIGMA1_h)*1000000)
disp(norm(Kfs*SIGMA1_h)*1000000)
disp(norm(Ksym*SIGMA1_h)*1000000)

SIGMA1_h    = reshape(SIGMA1_h, 6, MeshData.numnp);
SIGMA1_d    = N*SH1;
SIGMA1_d    = reshape(SIGMA1_d, 6, GridData.ndv);

titlestr    = { ...
    '\Sigma_{xx}^h'; ...
    '\Sigma_{yy}^h'; ...
    '\Sigma_{zz}^h'; ...
    '\Sigma_{yz}^h'; ...
    '\Sigma_{xz}^h'; ...
    '\Sigma_{xy}^h'};
figure(1)
for i = 1:1:6
    subplot(2,3,i)
    plot(SIGMA1_h(i,1:3)*1000000, 'b^')
    title(titlestr{i,1})
    xlabel('node number')
    ylabel('stress (MPa)')
end

titlestr    = { ...
    '\Sigma_{xx}^d'; ...
    '\Sigma_{yy}^d'; ...
    '\Sigma_{zz}^d'; ...
    '\Sigma_{yz}^d'; ...
    '\Sigma_{xz}^d'; ...
    '\Sigma_{xy}^d'};
figure(10)
for i = 1:1:6
    subplot(2,3,i)
    plot(SIGMA1_d(i,1:3)*1000000, 'b^')
    title(titlestr{i,1})
    xlabel('dv number')
    ylabel('stress (MPa)')
end

%%% METHOD 2 - SOLVE IN STEPS
% FIRST SOLVE THE MICROSCALE
SH2 = lsqlin(A, Fvec, ...
    [], [], ... 
    [], [], ...
    LB, UB, ...
    IG);

% SOLVE MACROSCALE USING THE MICROSCALE RESULTS
RHS         = Q*N*SH2;
SIGMA2_h    = K\RHS;

% CONCAT SOLUTION
SOL2    = [SIGMA2_h; SH2];

SIGMA2_h    = SOL2(1:nK,1);
SH2         = SOL2(nK+1:end,:);

disp(norm(Keq*SIGMA2_h)*1000000)
disp(norm(Kfs*SIGMA2_h)*1000000)
disp(norm(Ksym*SIGMA2_h)*1000000)

SIGMA2_h    = reshape(SIGMA2_h, 6, MeshData.numnp);
SIGMA2_d    = N*SH2;
SIGMA2_d    = reshape(SIGMA2_d, 6, GridData.ndv);

titlestr    = { ...
    '\Sigma_{xx}^h'; ...
    '\Sigma_{yy}^h'; ...
    '\Sigma_{zz}^h'; ...
    '\Sigma_{yz}^h'; ...
    '\Sigma_{xz}^h'; ...
    '\Sigma_{xy}^h'};
figure(2)
for i = 1:1:6
    subplot(2,3,i)
    plot(SIGMA2_h(i,1:3)*1000000, 'b^')
    title(titlestr{i,1})
    xlabel('node number')
    ylabel('stress (MPa)')
end

titlestr    = { ...
    '\Sigma_{xx}^d'; ...
    '\Sigma_{yy}^d'; ...
    '\Sigma_{zz}^d'; ...
    '\Sigma_{yz}^d'; ...
    '\Sigma_{xz}^d'; ...
    '\Sigma_{xy}^d'};
figure(20)
for i = 1:1:6
    subplot(2,3,i)
    plot(SIGMA2_d(i,1:3)*1000000, 'b^')
    title(titlestr{i,1})
    xlabel('dv number')
    ylabel('stress (MPa)')
end

%%% METHOD 3 - JUST SOLVE MICRO WITH APPLIED STRESS
Fmacro3 = repmat(InputData.SIMGA_APPLIED, 1, GridData.ndv);
Fmacro3 = Fmacro3(:);

F       = [Fmacro3; Fvec];

LHS = [N; A];

SH3 = lsqlin(LHS, F, ...
    [], [], ... 
    [], [], ...
    LB, UB, ...
    IG);

disp(norm(Keq*Q*N*SH3)*1000000)
disp(norm(Kfs*Q*N*SH3)*1000000)
disp(norm(Ksym*Q*N*SH3)*1000000)

SIGMA3_d    = N*SH3;
SIGMA3_d    = reshape(SIGMA3_d, 6, GridData.ndv);

titlestr    = { ...
    '\Sigma_{xx}^d'; ...
    '\Sigma_{yy}^d'; ...
    '\Sigma_{zz}^d'; ...
    '\Sigma_{yz}^d'; ...
    '\Sigma_{xz}^d'; ...
    '\Sigma_{xy}^d'};
figure(30)
for i = 1:1:6
    subplot(2,3,i)
    plot(SIGMA3_d(i,1:3)*1000000, 'b^')
    title(titlestr{i,1})
    xlabel('dv number')
    ylabel('stress (MPa)')
end

disp('Saving solution ...')
disp(FNAME_SOL)
save(PFNAME_SOL, 'SOL2')