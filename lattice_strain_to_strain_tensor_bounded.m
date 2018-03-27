clear all
close all
clc

% ARB STRAIN STATE 
data_file       = 'SampleListOfMeasuredStrains.xlsx';
[NUM,TXT,RAW]   = xlsread(data_file);

% EXPERIMENT
ome = [  0; 15; 30; 45; -45; 0;  0]; % ROTATION ABOUT YL
chi = [  0;  0;  0;  0;  0;  8; -8]; % ROTATION ABOUT XL

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
    0 sind(-th(1)) cosd(-th(1))]*y;

% SCATTERING VECTOR FOR HORZ DET
q(:,2)  = [ ...
    cosd(th(2)) 0 sind(th(2)); ...
    0 1 0; ...
    -sind(th(2)) 0 cosd(th(2))]*x;

% PROJECTION VECTOR FOR qv and qh
g(:,1)  = ProjectionVector(q(:,1)); %% V
g(:,2)  = ProjectionVector(q(:,2)); %% H

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
    
    %%% ARRANGE PROJECTION OPERATOR AND MEASURED STRAIN
    T   = VectorizedCOBMatrix(R);
    P   = [P; ...
        g(:,2)'*T; ...
        g(:,1)'*T; ...
        ];
end

A   = P'*P;
disp(sprintf('condition number for A : %3.1f', cond(A)))
disp(sprintf('condition number for P : %3.1f', cond(P)))

for iii = 1:1:size(NUM)
    disp(iii)
    xyz_sam(iii,:)  = NUM(iii,1:3);
    meas_strain_vec = NUM(iii,4:end)';
    
    lb  = [-inf -inf -inf -1e-3 -inf -inf]';
    ub  = [inf inf inf 1e-3 inf inf]';
    
    %%%
    % SQUARIZED LINEAR SYSTEM OF EQUATION THEN SOLVE
    b   = P'*meas_strain_vec;
    
    %%%
    % SOLVE BY LINEAR SOLVER1
%     [strain_vec(:,iii), flag(iii), relres(iii)] = lsqr(A, b);
%     resn(iii)   = norm(b - A*strain_vec);
    
    %%%
    % SOLVE BY LINEAR SOLVER2
    [strain_vec(:,iii), resn(iii), ~, flag(iii)]    = lsqlin(A, b, [], [], [], [], lb, ub);
end

% disp('inversion of itself+noise')
% MatrixOfStressStrainVectorInVM(strain_vecn)
figure(1)
for i = 1:1:6
    subplot(2,3,i)
    if i >= 4
        scatter(xyz_sam(:,2), strain_vec(i,:)./sqrt(2), 20, resn, 'filled')
    else
        scatter(xyz_sam(:,2), strain_vec(i,:), 20, resn, 'filled')
    end
    axis([min(xyz_sam(:,2)) max(xyz_sam(:,2)) -6e-3 6e-3])
end

