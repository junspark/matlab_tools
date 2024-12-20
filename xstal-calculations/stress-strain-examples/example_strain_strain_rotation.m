clear all
close all
clc

% addpath(genpath('/home/beams/PARKJS/matlab/matlab_tools'));
addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USING '11-22-33-23-13-12' ORDER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROTATOIN TAKES CRYSTAL FRAME TO SAMPLE FRAME
% R*c = s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rod     = [0.1657, 0.0429, 0.1036]';
quat    = QuatOfRod(rod);
R       = RMatOfQuat(quat);
T       = VectorizedCOBMatrix(R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARB STRAIN STATE IN XSTAL FRAME
e11 = 0.11;
e12 = -0.12;    % THIS IS EPSILON / NOT GAMMA
e13 = -0.13;    % THIS IS EPSILON / NOT GAMMA
e22 = 0.25;
e33 = -0.31;
e23 = 0.23;     % THIS IS EPSILON / NOT GAMMA

% BUILD STRAIN
strain_mtx  = [ ...
    e11 e12 e13; ...
    e12 e22 e23; ...
    e13 e23 e33; ...
    ];

strain_vec_in_c = VectorOfStressStrainMatrixInVM(strain_mtx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRYSTAL ELASTICITY
% SHEAR RELATES ij-COMPONENT OF STRESS TO 2*epsilon_ij (gamma_ij) COMPONENT OF STRAIN 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE 1
% CALCULATE STRESS IN CRYSTAL FRAME THEN TAKE STRESS TO SAMPLE FRAME
stress_vec_in_c     = C*strain_vec_in_c;
stress_vec_in_s1    = T*stress_vec_in_c;
stress_mtx_in_s1    = MatrixOfStressStrainVectorInVM(stress_vec_in_s1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CASE 2
% TAKE STRAIN TO SAMPLE FRAME THEN CALCULATE STRESS IN SAMPLE FRAME
strain_vec_in_s     = T*strain_vec_in_c;
stress_vec_in_s2    = T*C*T'*strain_vec_in_s;
stress_mtx_in_s2    = MatrixOfStressStrainVectorInVM(stress_vec_in_s2)