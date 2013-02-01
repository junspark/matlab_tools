clear all
close all
clc

%%% SOLVER COMBINING K Q A AND N
% RUN ONCE AND SAVE RESULT OR ELSE WASTE HUGE TIME
% MIGHT HAVE TO MIGRATE TO BIGGER SOLVER
load('/home/jspark/Documents/work/prjResidualStress/3D.MultiscaleCode2/A1234.K.Q.MLS.100,0.1,0.1,7.mat')
load('/home/jspark/Documents/work/Data/APS201103/GE/A1234.A.mat')
load('/home/jspark/Documents/work/Data/APS201103/GE/A1234.N.mat')
load('/home/jspark/Documents/work/Data/APS201103/GE/A1234.Fvec.mat')
load('/home/jspark/Documents/work/Data/APS201103/GE/A1234.INVIG.mat')
[nA, mA]    = size(A);
[nK, mK]    = size(K);
[nN, mN]    = size(N);
[nF, mF]    = size(Fvec);
[nIG, mIG]  = size(IG);

%%% FIRST SOLVE THE MICROSCALE
W   = A\Fvec;
W   = lsqlin(A, Fvec, ...
    [], [], ... 
    [], [], ...
    LB, UB, ...
    W);
% norm(Fvec - A*W)

%%% SOLVE MACROSCALE USING
%%% TAKING THE MICROSCALE RESULTS 
RHS     = Q*N*W;
SIGMA   = K\RHS;

clear K Q N A RHS

SH  = [SIGMA; W];

% Fvec_calc   = LHS*SH;
clear LHS Fvec 
save('A1234.SH.NEWRESULTS.mat', 'SH')
disp('done')