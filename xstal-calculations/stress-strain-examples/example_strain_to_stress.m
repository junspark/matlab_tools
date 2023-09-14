clear all
close all
clc

% addpath(genpath('/home/beams/PARKJS/matlab/matlab_tools'));
addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE THAT STRAIN IS ALREADY IN THE CRYSTAL FRAME
% USING '11-22-33-23-13-12' ORDER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARBITRARY STRAIN STATE
e11 = 0.11;
e12 = -0.12;    % THIS IS EPSILON / NOT GAMMA
e13 = -0.13;    % THIS IS EPSILON / NOT GAMMA
e22 = 0.25;
e33 = -0.31;
e23 = 0.23;     % THIS IS EPSILON / NOT GAMMA

% ELASTICITY
%%% CUBIC
% SHEAR RELATES ij-COMPONENT OF STRESS TO 2*epsilon_ij (gamma_ij) COMPONENT
% OF STRAIN (~Copper)
c11 = 168;
c12 = 121;
c44 = 75;

c(1)    = c11;
c(2)    = c12;
c(3)    = c44;

C = BuildElasticityMatrix(c, 'Symmetry', 'cubic')

C_conventional  = [ ...
    c11 c12 c12 0 0 0 ;...
    c12 c11 c12 0 0 0 ; ...
    c12 c12 c11 0 0 0 ; ...
    0 0 0 c44 0 0 ; ...
    0 0 0 0 c44 0 ; ...
    0 0 0 0 0 c44]

%%% STRAIN MATRIX IN EPSILON
strain_mtx  = [ ...
    e11 e12 e13; ...
    e12 e22 e23; ...
    e13 e23 e33; ...
    ];

%%% SHEAR STRAIN WITH sqrt(2) FACTOR
strain_vec  = VectorOfStressStrainMatrixInVM(strain_mtx);
%%% SHEAR STRAIN IN GAMMA
strain_conventional = [e11 e22 e33 2*e23 2*e13 2*e12]';           

%%% THIS IS STRESS WITH sqrt(2) FACTOR IN THE SHEARS
stress_vec  = C*strain_vec;

%%% THIS IS STRESS MATRIX IN CONVENTIONAL NOTATION (SAME AS BELOW)
stress_mtx  = MatrixOfStressStrainVectorInVM(stress_vec);

%%% THIS IS THE USUAL SIGMA = C*EPSILON
stress_conventional = C_conventional*strain_conventional;

%%% HEXAGONAL
% SHEAR RELATES ij-COMPONENT OF STRESS TO 2*epsilon_ij (gamma_ij) COMPONENT
% OF STRAIN (~Titanium)
c11 = 162;
c33 = 181;
c44 = 45;
c12 = 92;
c13 = 69;

c(1)    = c11;
c(2)    = c33;
c(3)    = c12;
c(4)    = c13;
c(5)    = c44;

C = BuildElasticityMatrix(c, 'Symmetry', 'hexagonal');

C_conventional  = [ ...
    c11 c12 c13 0 0 0 ;...
    c12 c11 c13 0 0 0 ; ...
    c13 c13 c33 0 0 0 ; ...
    0 0 0 (c11 - c12)/2 0 0 ; ...
    0 0 0 0 c44 0 ; ...
    0 0 0 0 0 c44];

% BUILD STRAIN
strain_mtx  = [ ...
    e11 e12 e13; ...
    e12 e22 e23; ...
    e13 e23 e33; ...
    ];
strain_vec  = VectorOfStressStrainMatrixInVM(strain_mtx);
strain_conventional = [e11 e22 e33 2*e23 2*e13 2*e12]';

stress_vec  = C*strain_vec;
stress_mtx  = MatrixOfStressStrainVectorInVM(stress_vec)

stress_conventional = C_conventional*strain_conventional