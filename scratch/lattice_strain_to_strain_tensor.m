clear all
close all
% clc

% ARB STRAIN STATE 1
e11 = 0.11;
e12 = -0.12;
e13 = -0.13;
e22 = 0.25;
e33 = -0.31;
e23 = 0.23;

% ARB STRAIN STATE 2
% e11 = 0.11;
% e12 = -0.12;
% e13 = 0;
% e22 = 0.25;
% e33 = 0;
% e23 = 0;

strain_mtx  = [ ...
    e11 e12 e13; ...
    e12 e22 e23; ...
    e13 e23 e33; ...
    ];
strain_vec  = VectorOfStressStrainMatrixInVM(strain_mtx);

% EXPERIMENT
% ome = [-45; 0; 15; 30; 45;  0; 0; 0]; % ROTATION ABOUT YL
% chi = [  0; 0;  0;  0;  0; -8; 0; 8]; % ROTATION ABOUT XL
% ome = [  0; 45; 90;  0;  0; 45]; ; % ROTATION ABOUT YL
% chi = [  0;  0;  0; 45; 90; 45]; % ROTATION ABOUT XL
% ome = [  0;   0;  0]; % ROTATION ABOUT YL
% chi = [  0;  45; 90]; % ROTATION ABOUT XL
ome = zeros(72,1); % ROTATION ABOUT YL
chi = 0:5:355; % ROTATION ABOUT XL (CAN BE ETA WITH AREA DETECTOR)

% VERT / HORZ DETECTOR
tth(1)  = 4.99995; % VERT
tth(2)  = 4.90353; % HORZ

th      = tth./2;
% th      = zeros(1,2);

x       = [1 0 0]';
y       = [0 1 0]';

% SCATTERING VECTOR FOR VERT DET
q(:,1)  = [ ...
    1 0 0; ...
    0 cosd(-th(1)) -sind(-th(1)); ...
    0 sind(-th(1)) cosd(-th(1))]*[0 1 0]';

% SCATTERING VECTOR FOR HORZ DET
q(:,2)  = [ ...
    cosd(th(2)) 0 sind(th(2)); ...
    0 1 0; ...
    -sind(th(2)) 0 cosd(th(2))]*[1 0 0]';

% PROJECTION VECTOR FOR qv and qh
g(:,1)  = ProjectionVector(q(:,1));
g(:,2)  = ProjectionVector(q(:,2));

P   = [];
for i = 1:1:length(ome)
    disp(sprintf('%d | %3.1f | %3.1f', i, ome(i), chi(i)))
    come    = cosd(ome(i));
    some    = sind(ome(i));
    
    cchi    = cosd(chi(i));
    schi    = sind(chi(i));
    %%% GOING FROM SAMPLE FRAME TO LAB FRAME
    Ry   = [ ...
        come 0 some; ...
        0 1 0; ...
        -some 0 come];
    
    %%%%%
    Rz   = [ ...
        cchi -schi 0; ...
        schi cchi 0; ...
        0 0 1];
    
    %%% NOTE THAT THIS IS SAYING 
    %%% ROTATION ABOUT Y THEN 
    %%% ROTATION ABOUT THE ROTATED Z
    R   = Rz*Ry;
    
    %%% 
    T   = VectorizedCOBMatrix(R);
    P   = [P; ...
        g(:,1)'*T; ...
        g(:,2)'*T; ...
        ];
    
    % R*strain_mtx*R'
    % MatrixOfStressStrainVectorInVM(T*strain_vec)
    % return
end

%%%
meas_strain = P*strain_vec;
a1  = -2e-4;
a2  = 2e-4;

noise   = a1+ (a2-a1).*rand(length(meas_strain),1);

% SQUARIZED LINEAR SYSTEM OF EQUATION THEN SOLVE
A   = P'*P;
b0  = P'*meas_strain;
bn  = P'*(meas_strain + noise);

strain_vec0	= A\b0;
strain_vecn = A\bn;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(sprintf('condition number for A : %3.1f', cond(A)))
disp(sprintf('condition number for P : %3.1f', cond(P)))
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

disp('answer')
strain_mtx
disp('inversion of itself')
MatrixOfStressStrainVectorInVM(strain_vec0)
disp('inversion of itself+noise')
MatrixOfStressStrainVectorInVM(strain_vecn)
