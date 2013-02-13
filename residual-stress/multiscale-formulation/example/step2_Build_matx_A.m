clear all
close all
clc

%%% INPUT DATA GENERATED IN STEP0
InputData   = load('example_rsinput.mat');

%%%
% FIBER IS A TUBE IN ORIENTATION SPACE DUE TO OMEGA AND ETA STEP SIZES
% COMPUTE WEIGHTS FOR FIBER TUBE
w = FiberTubeWeight(InputData.deta, InputData.domega);

ct  = 0;
for i = 1:1:InputData.ndv
    %%%
    % IN THIS CASE, THE SCATTERING VECTORS ARE ALWAYS THE SAME
    % SO WE USE THE A0 FROM STEP1 AS IS.
    % IN THE CASE OF CONICAL, WE NEED TO TAKE OUT 
    % ROWS FROM A0 GENERATED IN STEP1.
    if ct == 6
        ct  = 0;
    end
    if InputData.UseLSDF
        disp('Generate LSDF <> LS ...')
        f   = ['BaseA.LSDF2LS.DV', num2str(ct), '.mat'];
        pf  = fullfile(InputData.PNAME_SOLUTION, f);
    else
        disp('Generate SODF <> LS ...')
        f   = ['BaseA.SODF2LS.DV', num2str(ct), '.mat'];
        pf  = fullfile(InputData.PNAME_SOLUTION, f);
    end
    
    if exist(pf,'file')
        disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
        disp(f)
        disp('exists')
    else
        disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
        disp(f)
        disp('does not exist')
    end
    ct  = ct + 1;
end