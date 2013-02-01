clear all
close all
clc

%%% INPUT FILES FROM HERE ARE USE UP TO STEP3

%%% NAME OF INPUT FILE
InputFileName   = 'inputdata.A1.mat';       %% BE SURE TO CHANGE LINE 18 AS WELL

%%% INVERSION OPTION
inputdata.nHarm         = 10;       % NUMBER OF SHM USED
inputdata.Use20xMesh    = 1;        % FLAG: USE PREGENERATED SHM FOR CONSISTENCY
inputdata.UsePM         = 1;        % FLAG: PREMULTIPLIER
inputdata.eqq_calc_fig  = 1;        % FLAG: CALCULATE eqq AFTER INVERSION

%%% PATH DESIGNATION
inputdata.PNAME_SAM         = '/home/jspark/Documents/work/Data/APS201103/GE'; % data set path
inputdata.PNAME_LS_FILTER   = fullfile(inputdata.PNAME_SAM, 'A1');
inputdata.PNAME_SOLUTION    = '/home/jspark/Documents/work/prjResidualStress/3D.MultiscaleCode3/SOLUTION/'; % solution path

%%% DV grid
inputdata.x_grid    = [0.000];
inputdata.y_grid    = [-0.500; -0.900; -1.450; -1.950; -2.400; -2.950];
inputdata.z_grid    = [-0.525; -0.375; -0.225; -0.075; +0.075; +0.225; +0.375; +0.525];

%%% ETA
inputdata.eta_bins  = 72;                               % number of azimuthal bins
inputdata.eta_ini   = -360/inputdata.eta_bins/2;        % start azimuth (min edge of bin) in degrees
inputdata.eta_fin   = 360-360/inputdata.eta_bins/2;     % stop  azimuth (max edge of bin) in degrees
inputdata.eta_grid  = 0:360/inputdata.eta_bins:inputdata.eta_fin;

%%% OMEGA
inputdata.omega_bins    = 25;           % number of azimuthal bins
inputdata.omega_ini     = -60;          % start azimuth (min edge of bin) in degrees
inputdata.omega_fin     = 60;           % stop  azimuth (max edge of bin) in degrees
inputdata.omega_grid    = linspace(inputdata.omega_ini,inputdata.omega_fin,inputdata.omega_bins);

%%% X-RAY ENERGY
inputdata.energy        = 61.999;       % keV
inputdata.wavelength    = keV2Angstrom(inputdata.energy);  % wavelength (Angstrom)

%%% PROJECTION OPERATOR PARAMETERS
inputdata.sym       = CubSymmetries;    % xstal symmetry
inputdata.div       = 1000;             % number of divisions per fiber
inputdata.useA0     = 1;                % use baseline A0 to generate A matrices (need to run step0a)
inputdata.saveA     = 1;                % flag: save projection matrix A per DV
inputdata.useLSDF   = 0;                % use baseline A0 to generate A matrices (need to run step0a)

%%% FILENAMES
inputdata.FNAME_A       = 'A1234.A.mat';
inputdata.FNAME_N       = 'A1234.N.mat';
inputdata.FNAME_L       = 'A1234.L.mat';
inputdata.FNAME_Nz      = 'A1234.Nz.mat';
inputdata.FNAME_Nxy     = 'A1234.Nxy.mat';
inputdata.FNAME_Fvec    = 'A1234.Fvec.mat';
inputdata.FNAME_LS      = 'A1234.LSALL.mat';
inputdata.FNAME_Fmk     = 'A1234.Fmk.mat';
inputdata.FNAME_Re      = 'A1234.Re.mat';
inputdata.FNAME_RANK    = 'A1234.RANK.mat';
inputdata.FNAME_INVIG   = 'A1234.INVIG.mat';
inputdata.FNAME_GRID    = 'A1234.DV.GRID.mat';

%%% HKL LATTICE STRAIN TO USE
inputdata.hkls  = [1 1 1; 2 0 0; 2 2 0];

%%% ODF
%%% NOTE THAT FRMESH ASSOCIATED WITH THE ODF IS USED TO REPRESENT SODF
inputdata.pname_odf = '/home/jspark/Documents/work/MaterialProperty/LSHR.Ni/Texture';
inputdata.fname_odf = 'texture.odf.mat';

%%%
% MATL PROPERTY
% SX-STIFFNESS
C11 = 245.25e3;     % MPa
C12 = 155.16e3;     % MPa
C44 = 125e3;        % MPa CAREFUL FOR FACTOR OF TWO!!! SIG_XY = C44 * GAMMA_XY
C   = [C11; C12; 2*C44];
S   = C2S(C);
inputdata.C = C;
inputdata.S = S;

inputdata.lattparms = 3.5924;  % LSHR

disp('generating input file')
save(InputFileName, 'inputdata')