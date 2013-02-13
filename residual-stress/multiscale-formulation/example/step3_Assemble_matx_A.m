clear all
close all
clc

%%% INPUT DATA GENERATED IN STEP0
InputData   = load('example_rsinput.mat');

%%% VARIOUS PFNAMES
PFNAME_A        = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_A);
PFNAME_N        = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_N);
PFNAME_L        = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_L);
PFNAME_Nz       = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_Nz);
PFNAME_Nxy      = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_Nxy);
PFNAME_Fvec     = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_Fvec);
PFNAME_LS       = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_LS);
PFNAME_Fmk      = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_Fmk);
PFNAME_Re       = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_Re);
PFNAME_RANK     = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_RANK);
PFNAME_INVIG    = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_INVIG);
PFNAME_Ksym     = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_Ksym);
PFNAME_Kfs      = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_Kfs);
PFNAME_Keq      = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_Keq);

%%%
% PREP FRMESH
frmesh  = InputData.frmesh;
nnode   = frmesh.numind;
RMats   = RMatOfQuat(QuatOfRod(frmesh.crd(:,1:nnode)));

%%% GENERATE SH MODES
if InputData.UseProjectedHarmonics
    disp('projecting standard harmonics to frmesh ...')
    load wshex12x
    dh  = zeros(nnode, InputData.nHarm);
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
DH	= zeros(6*nnode,6*InputData.nHarm);

DH(1:6:end-5,1:InputData.nHarm)                     = dh;           % 11
DH(2:6:end-4,InputData.nHarm+1:2*InputData.nHarm)   = dh;           % 22
DH(3:6:end-3,2*InputData.nHarm+1:3*InputData.nHarm) = dh;           % 33
DH(4:6:end-2,3*InputData.nHarm+1:4*InputData.nHarm) = dh;           % 12
DH(5:6:end-1,4*InputData.nHarm+1:5*InputData.nHarm) = dh;           % 13
DH(6:6:end,5*InputData.nHarm+1:6*InputData.nHarm)   = dh;           % 23

DH	= sparse(DH);

%%%
% SX-STIFFNESS
C   = InputData.C;
S   = InputData.S;

%%% BUILD BigC or BigS - STIFFNESS / COMPLIANCE MATRIX
if InputData.UseLSDF
    disp('generating BigC ...')
    BigC    = zeros(nnode*6,nnode*6);
    for p = 1:1:nnode
        T   = VectorizedCOBMatrix(RMats(:,:,p));
        ri  = 6*(p-1) + 1;
        rf  = 6*p;
        ci  = 6*(p-1) + 1;
        cf  = 6*p;
        BigC(ri:rf,ci:cf)   = T*C*T';
    end
    BigC    = sparse(BigC);
else
    disp('Using SODF - compliance was included in A0 matrix ...')
end

%%%
% BUILD WEIGHTING MATRIX
BigL2IP = zeros(nnode*6,nnode*6);
for p = 1:1:nnode
    st  = 6*(p-1)+1;
    ed  = 6*p;
    
    BigL2IP(st,1:6:end-5)   = frmesh.l2ip(p,:);
    BigL2IP(st+1,2:6:end-4) = frmesh.l2ip(p,:);
    BigL2IP(st+2,3:6:end-3) = frmesh.l2ip(p,:);
    BigL2IP(st+3,4:6:end-2) = frmesh.l2ip(p,:);
    BigL2IP(st+4,5:6:end-1) = frmesh.l2ip(p,:);
    BigL2IP(st+5,6:6:end-0) = frmesh.l2ip(p,:);
end
BigL2IP = sparse(BigL2IP);

%%%
% BUILD A
%%%

%%% INITIALIZE VARIABLES
q   = [];
eqq = [];
eqq_mk  = [];

BigA        = [];
BigNMat     = [];
BigNMatz    = [];
BigNMatij   = [];
BigL        = [];
Fvec        = [];

IG  = [];
UB  = [];
LB  = [];

q_last  = 0;
HarmNum = 0;
SijNum  = 0;

% STRESSES
SIGMA11 = zeros(InputData.ndv,1);
SIGMA22 = zeros(InputData.ndv,1);
SIGMA33 = zeros(InputData.ndv,1);
SIGMA12 = zeros(InputData.ndv,1);
SIGMA13 = zeros(InputData.ndv,1);
SIGMA23 = zeros(InputData.ndv,1);
STRESS  = [];

% WEIGHTED ERROR
Re      = cell(size(InputData.hkls,1)+1,1);
% RANK
Rk      = zeros(InputData.ndv,4);

tic
PosNum  = 0;
for i = 1:1:InputData.ndv
    if PosNum == 6
        PosNum  = 0;
    end
    %%% LOAD ODF
    disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
    disp('loading odf for dv ...')
    odf_data    = load(InputData.PFNAME_ODF{i,1});
    % PlotFR(frmesh, odf_data, 'Symmetries', 'hexagonal')
    
    sum_l2  = sum(frmesh.l2ip*odf_data);
    
    %%%
    % ASSEMBLE BigODF MATRIX
    % NEEDED FOR SUMMATION OPERATOR N TO SUM STRESSES OVER OS
    BigODF              = zeros(6*nnode,6);
    BigODF(1:6:end-5,1) = odf_data;
    BigODF(2:6:end-4,2) = odf_data;
    BigODF(3:6:end-3,3) = odf_data;
    BigODF(4:6:end-2,4) = odf_data;
    BigODF(5:6:end-1,5) = odf_data;
    BigODF(6:6:end,6)   = odf_data;
    BigODF              = sparse(BigODF);
    
    %%% LOAD SCATTERING VECTORS AND LATTICE STRAINS PER DV
    disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
    q_DV    = [];
    eqq_DV  = [];
    q_size_DV   = [];
    q_last_DV   = 0;
    for j = 1:1:InputData.nhkls
        fname_LS_DATA   = [ ...
            InputData.FNAME_LS_DATA_ROOT, ...
            num2str(InputData.hkls(j,1)), ...
            num2str(InputData.hkls(j,2)), ...
            num2str(InputData.hkls(j,3)), ...
            '.data'];
        disp(['dv num : ', num2str(PosNum), ' - loading ', fname_LS_DATA])
        
        pfname_LS_DATA  = fullfile(InputData.PNAME_LS_DATA{i,1}, fname_LS_DATA);
        
        LS_data	= load(pfname_LS_DATA);
        hkl     = ones(size(LS_data,1),1)*InputData.hkls(j,:);
        
        q_DV    = [q_DV; LS_data(:,1:3)];
        eqq_DV  = [eqq_DV; LS_data(:,4)];
        
        idx_q_DV    = size(LS_data,1) + q_last_DV;
        
        q_size_DV   = [q_size_DV; ...
            i InputData.hkls(j,:) idx_q_DV];
        
        q_last_DV   = q_size_DV(end,5);
        
        % FOR ALL DVs
        q   = [q; LS_data(:,1:3)];
        eqq = [eqq; LS_data(:,4)]; % RHS OF THE SOLVER
        
        idx_q   = size(LS_data,1) + q_last;
        eqq_mk  = [eqq_mk; ...
            i InputData.hkls(j,:) idx_q];
        
        q_last  = eqq_mk(end,5);
    end
    
    %%%
    % IN THIS CASE, THE SCATTERING VECTORS ARE ALWAYS THE SAME
    % SO WE USE THE A0 FROM STEP1 AS IS.
    % IN THE CASE OF CONICAL, WE NEED TO TAKE OUT 
    % ROWS FROM A0 GENERATED IN STEP1.
    if InputData.UseLSDF
        disp('Generate LSDF <> LS ...')
        f   = ['BaseA.LSDF2LS.DV', num2str(PosNum), '.mat'];
        pf  = fullfile(InputData.PNAME_SOLUTION, f);
    else
        disp('Generate SODF <> LS ...')
        f   = ['BaseA.SODF2LS.DV', num2str(PosNum), '.mat'];
        pf  = fullfile(InputData.PNAME_SOLUTION, f);
    end
    
    if exist(pf,'file') == 2
        disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
        disp(f)
        disp('exists')
        load(pf)
    else
        disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
        disp(f)
        disp('does not exist')
        A0  = sparse(size(q_DV,1),6*nnode);
    end
    
    %%%
    % ASSEMBLE A INTO BigA
    A0  = A0*DH;
    if InputData.UsePreMultiplier
        ATA     = A0'*A0;
        
        ri  = size(BigA,1) + 1;
        rf  = size(BigA,1) + InputData.nHarm*6;
        ci  = InputData.nHarm*6*(i-1) + 1;
        cf  = InputData.nHarm*6*i;
        BigA(ri:rf, ci:cf)  = ATA;
        
        ri  = size(BigL,1) + 1;
        rf  = size(BigL,1) + size(A0,1);
        ci  = InputData.nHarm*6*(i-1) + 1;
        cf  = InputData.nHarm*6*i;
        BigL(ri:rf, ci:cf)  = A0;
        
        Fvec    = [Fvec; A0'*eqq_DV];
    else
        ri  = size(BigA,1) + 1;
        rf  = size(BigA,1) + size(A0,1);
        ci  = InputData.nHarm*6*(i-1) + 1;
        cf  = InputData.nHarm*6*i;
        BigA(ri:rf, ci:cf)  = A0;
        
        Fvec    = [Fvec; eqq_DV];
    end
    
    %%%
    % SUMMATION OPERATOR TO SUM STRESS OVER ORIENTATION
    % SHEARS STILL HAVE SQRT(2)
    if InputData.UseLSDF
        Nmat	= BigODF'*BigL2IP*BigC*DH/sum_l2;
    else
        Nmat	= BigODF'*BigL2IP*DH/sum_l2;
    end
    
    %%% GENERATE SIGMA_APPLIED TO BE USED FOR IG GENERATION
    %%% USE SIMPLE INVERSION TO GENERATE APPRAENT MACRO STRESS
    if InputData.UsePreMultiplier
        X_DV        = ATA\(A0'*eqq_DV);     %%% SIMPLEST INV
    else
        X_DV        = A0\eqq_DV;            %%% SIMPLEST INV
    end
    
    % SHEARS STILL HAVE SQRT(2)
    % SIGMA_APPLIED   = Nmat*X_DV;
    SIGMA_APPLIED   = InputData.SIMGA_APPLIED;
    
    %%%
    % PERFORM SIMPLE INVERSION
    %%%
    %%% FIRST GENERATE IG, UB, LB
    %%% CALCULATE INV-IG
    S_MACRO         = AggregateElasticityMatrix(frmesh, odf_data, S);
    EPSILON_APPLIED = S_MACRO*SIGMA_APPLIED;
    clear ig ub lb
    [ig, lb, ub]    = CalculateSPFInversionIG(frmesh, dh, ...
        EPSILON_APPLIED, C);

%     k2  = 0.5;
%     ub  = ig + abs(ig)*k2;
%     lb  = ig - abs(ig)*k2;
%     
%     iub = abs(ig) < 0.5;
%     ilb = abs(ig) < 0.5;
%     
%     ub(iub) = +1.0;
%     lb(ilb) = -1.0;

    IG  = [IG; ig];
    UB  = [UB; ub];
    LB  = [LB; lb];
    
    %%% PLOT IG
    SHIG    = reshape(ig', InputData.nHarm, 6);
    SHIG_11 = SHIG(:,1);
    SHIG_22 = SHIG(:,2);
    SHIG_33 = SHIG(:,3);
    SHIG_12 = SHIG(:,4);
    SHIG_13 = SHIG(:,5);
    SHIG_23 = SHIG(:,6);
    
    SODF11_IG   = dh*SHIG_11;
    SODF22_IG   = dh*SHIG_22;
    SODF33_IG   = dh*SHIG_33;
    SODF12_IG   = dh*SHIG_12./sqrt(2);
    SODF13_IG   = dh*SHIG_13./sqrt(2);
    SODF23_IG   = dh*SHIG_23./sqrt(2);
    
    % PlotFR(frmesh, SODF11_IG)
    % PlotFR(frmesh, SODF22_IG)
    % PlotFR(frmesh, SODF33_IG)
    % PlotFR(frmesh, SODF12_IG)
    % PlotFR(frmesh, SODF13_IG)
    % PlotFR(frmesh, SODF23_IG)
    
    %%% INVERT
    if InputData.UsePreMultiplier
        % X_DV        = ATA\(A0'*eqq_DV);    %%% SIMPLEST INV
        X_DV    = lsqlin(ATA, A0'*eqq_DV, ...
            [], [], ...
            [], [], ...
            lb, ub, ...
            ig);
        eqq_calc_DV = A0*X_DV;
    else
        % X_DV        = A0\eqq_DV;           %%% SIMPLEST INV
        X_DV    = lsqlin(A0, eqq_DV, ...
            [], [], ...
            [], [], ...
            lb, ub, ...
            ig);
        eqq_calc_DV = A0*X_DV;
    end
    
    %%% PLOT INV RESULT
    SODF    = DH*X_DV;
    SODF11  = SODF(1:6:end-5);
    SODF22  = SODF(2:6:end-4);
    SODF33  = SODF(3:6:end-3);
    SODF12  = SODF(4:6:end-2)./sqrt(2);
    SODF13  = SODF(5:6:end-1)./sqrt(2);
    SODF23  = SODF(6:6:end-0)./sqrt(2);
    
    %                 PlotFR(frmesh, SODF11)
    %                 PlotFR(frmesh, SODF22)
    %                 PlotFR(frmesh, SODF33)
    %                 PlotFR(frmesh, SODF12)
    %                 PlotFR(frmesh, SODF13)
    %                 PlotFR(frmesh, SODF23)
    
    SIGMA       = Nmat*X_DV;        % UNITS DEPEND ON C INPUT
    SIGMA(4:6)  = SIGMA(4:6)./sqrt(2);
    
    MEAN_LS_DIFF(i,1)   = mean(abs(eqq_calc_DV - eqq_DV))/max(abs(eqq_DV));
    MAX_LS_DIFF(i,1)    = max(abs(eqq_calc_DV - eqq_DV))/max(abs(eqq_DV));
    
    disp('dv Sxx Syy Szz Syz Sxz Sxy FOR IG')
    disp('(UNITS IN MPa)')
    disp([PosNum round(SIGMA_APPLIED'.*1000000)])   %%% *1000000 TPa <> MPa
    disp('dv Sxx Syy Szz Syz Sxz Sxy SOL')
    disp('(UNITS IN MPa)')
    disp([PosNum round(SIGMA'.*1000000)])           %%% *1000000 TPa <> MPa
    disp('% mean LS diff')
    disp(mean(abs(eqq_calc_DV - eqq_DV)))
    disp('% max LS diff')
    disp(max(abs(eqq_calc_DV - eqq_DV)))
    
    SIGMA11(i,1)	= SIGMA(1,1);
    SIGMA22(i,1)	= SIGMA(2,1);
    SIGMA33(i,1)	= SIGMA(3,1);
    SIGMA12(i,1)	= SIGMA(6,1);
    SIGMA13(i,1)	= SIGMA(5,1);
    SIGMA23(i,1)	= SIGMA(4,1);
    STRESS(i,:)     = SIGMA';
    
    %%% 3D MACRO CODE TAKES SHEARS W/O sqrt(2)
    %%% 3D MACRO CODE STRESS ORDER 11 22 33 23 13 12
    Nmat([4 5 6],:) = Nmat([6 5 4],:)./sqrt(2);
    
    ri  = SijNum+1;
    rf  = SijNum+6;
    ci  = HarmNum+1;
    cf  = HarmNum+6*InputData.nHarm;
    BigNMat(ri:rf,ci:cf)    = Nmat;
    
    ri  = SijNum/2+1;
    rf  = SijNum/2+3;
    ci  = HarmNum+1;
    cf  = HarmNum+6*InputData.nHarm;
    BigNMatz(ri:rf,ci:cf)   = Nmat([3 5 6],:);
    BigNMatij(ri:rf,ci:cf)  = Nmat([1 2 4],:);
    
    %%%
    % MARGIN OF ERROR CALCULATION
    % SPF IMAGE GENERATION
    Re{end,1}   = [ ...
        Re{end,1}; ...
        InputData.dvx(i), ...
        InputData.dvy(i), ...
        InputData.dvz(i), ...
        ErrorRe(eqq_DV, eqq_calc_DV, 'Threshold', 0.0001)];
    Rk(i,:)     = [ ...
        InputData.dvx(i), ...
        InputData.dvy(i), ...
        InputData.dvz(i), ...
        rank(full(A0))/size(A0,2)];
    
    for j = 1:1:size(InputData.hkls,1)
        H   = InputData.hkls(j,1);
        K   = InputData.hkls(j,2);
        L   = InputData.hkls(j,3);
        
        if j == 1
            ri  = 1;
            rf  = q_size_DV(j,5);
        else
            ri  = q_size_DV(j-1,5) + 1;
            rf  = q_size_DV(j,5);
        end
        recalc_data = [q_DV(ri:rf,:) eqq_calc_DV(ri:rf)];
        org_data    = [q_DV(ri:rf,:) eqq_DV(ri:rf)];
        
        % CALCULATE ERROR MEASURES
        Re{j,1} = [Re{j,1}; ...
            InputData.dvx(i), ...
            InputData.dvy(i), ...
            InputData.dvz(i), ...
            ErrorRe(eqq_DV(ri:rf), eqq_calc_DV(ri:rf),'Threshold', 0.0001)];
        
        % SAVE RECALCULATED DATA
        if InputData.Plot_recalc_SPF || InputData.Save_recalc_SPF
            % SAVE RECALCULATED DATA
            if InputData.Save_recalc_SPF
                if InputData.UseLSDF
                    f   = [ ...
                        InputData.FNAME_LS_DATA_ROOT, ...
                        num2str(H), num2str(K), num2str(L), ...
                        '.LS.recalc.LSDF2LS.data' ...
                        ];
                else
                    f   = [ ...
                        InputData.FNAME_LS_DATA_ROOT, ...
                        num2str(H), num2str(K), num2str(L), ...
                        '.LS.recalc.SODF2LS.data' ...
                        ];
                end
                pf  = fullfile(InputData.PNAME_LS_DATA{i,1}, f);
                
                disp(['dv num : ', num2str(PosNum), ' - saving ', f])
                fid = fopen(pf, 'w');
                fprintf(fid, '%% RECALCULATED SPF DATA VIA INVERSION W/O MACRO BCS\n');
                fprintf(fid, '%% q eqq\n');
                fprintf(fid, '%2.6f %2.6f %2.6f %2.6f \n', recalc_data');
                fclose(fid);
            else
                disp('not saving recalculated strain data ...')
            end
            
            % SAVE RECALCULATED DATA SPF
            if InputData.Plot_recalc_SPF
                if InputData.UseLSDF
                    f   = [ ...
                        InputData.FNAME_LS_DATA_ROOT, ...
                        num2str(H), num2str(K), num2str(L), ...
                        '.LS.recalc.LSDF2LS' ...
                        ];
                else
                    f   = [ ...
                        InputData.FNAME_LS_DATA_ROOT, ...
                        num2str(H), num2str(K), num2str(L), ...
                        '.LS.recalc.SODF2LS' ...
                        ];
                end
                pf  = fullfile(InputData.PNAME_LS_DATA{i,1}, f);
                
                fig1    = figure(1000);
                PlotSPF(recalc_data(:,1:3), recalc_data(:,4), ...
                    'Title', 'SPF_{CALC}', ...
                    'DataRange', [-0.0056 0.0074])
                hold off
                
                disp(['dv num : ', num2str(PosNum), ' - saving ', f, '.png'])
                saveas(fig1, [pf, '.png'], 'png')
                disp(['dv num : ', num2str(PosNum), ' - saving ', f, '.fig'])
                saveas(fig1, [pf, '.fig'], 'fig')
                
                fig2    = figure(2000);
                PlotSPF(recalc_data(:,1:3), org_data(:,4) - recalc_data(:,4), ...
                    'Title', 'SPF_{EXP} - SPF_{CALC}')
                hold off
                
                disp(['dv num : ', num2str(PosNum), ' - saving ', f, '.png'])
                saveas(fig2, [pf, '.diff.png'], 'png')
                disp(['dv num : ', num2str(PosNum), ' - saving ', f, '.fig'])
                saveas(fig2, [pf, '.diff.fig'], 'fig')
                pause(0.25)
            end
        else
            disp('not saving recalculated strain pole figures ...')
        end
    end
    % pause
    close all
    
    % UPDATE COUNTERS
    HarmNum = HarmNum + 6*InputData.nHarm;
    SijNum  = SijNum + 6;
    PosNum  = PosNum + 1;
    
    % MEMORY MANAGEMENT
    clear A0 Nmat
    BigA        = sparse(BigA);
    BigNMat     = sparse(BigNMat);
    BigNMatz    = sparse(BigNMatz);
    BigNMatij   = sparse(BigNMatij);
    if InputData.UsePreMultiplier
        BigL    = sparse(BigL);
    end
end
toc

%%%
% ASSIGN VARIABLES AND SAVE FOR NEXT STEP
disp('generating DV inversion files')
A   = BigA;
N   = BigNMat;
Nz  = BigNMatz;
Nxy = BigNMatij;
L   = BigL;

clear BigA BigNMat BigNMatz BigNMatij BigL

save(PFNAME_A, 'A')
save(PFNAME_N, 'N')
save(PFNAME_L, 'L')
save(PFNAME_Nz, 'Nz')
save(PFNAME_Nxy, 'Nxy')
save(PFNAME_LS, 'eqq', 'q', 'eqq_mk')
save(PFNAME_Re, 'Re')
save(PFNAME_RANK, 'Rk')
save(PFNAME_INVIG , 'IG', 'UB', 'LB')
save(PFNAME_Fvec, 'Fvec')
return
%%%
%%% NEEDS TO BE UPDATED FOR THIS PROBLEM
% PERFORM INVERSION OF THE SYSTEM
X           = A\Fvec;
X_DV    = lsqlin(A, Fvec, ...
    [], [], ...
    [], [], ...
    LB, UB, ...
    IG);

%%%%
% PLOT SOME RESULTS TO CHECK

S11 = [];
S22 = [];
S33 = [];
S23 = [];
S13 = [];
S12 = [];
for i = 1:1:ny
    S11 = [S11; SIGMA11(1,i,:)];
    S22 = [S22; SIGMA22(1,i,:)];
    S33 = [S33; SIGMA33(1,i,:)];
    S12 = [S12; SIGMA12(1,i,:)];
    S13 = [S13; SIGMA13(1,i,:)];
    S23 = [S23; SIGMA23(1,i,:)];
end
max_S   = max([S11(:); S22(:); S33(:); S12(:); S13(:); S23(:)]);
min_S   = min([S11(:); S22(:); S33(:); S12(:); S13(:); S23(:)]);

SEFF    = sqrt(((S11 - S22).^2 + (S22 - S33).^2 + (S33 - S11).^2 + ...
    6*S12.*S12 + 6*S13.*S13 + 6*S23.*S23)./2);

figure(1)
surf(InputData.y_grid, InputData.z_grid, squeeze(S11))
title('\Sigma_{11}')
view([0 90])
caxis([min_S 700])
axis equal tight

figure(2)
surf(InputData.y_grid, InputData.z_grid, S22)
title('\Sigma_{22}')
view([0 90])
caxis([min_S 700])
axis equal tight

figure(3)
surf(InputData.y_grid, InputData.z_grid, S33)
title('\Sigma_{33}')
view([0 90])
caxis([min_S 700])
axis equal tight

figure(4)
surf(InputData.y_grid, InputData.z_grid, S12)
title('\Sigma_{12}')
view([0 90])
caxis([min_S 700])
axis equal tight

figure(5)
surf(InputData.y_grid, InputData.z_grid, S13)
title('\Sigma_{13}')
view([0 90])
caxis([min_S 700])
axis equal tight

figure(6)
surf(InputData.y_grid, InputData.z_grid, S23)
title('\Sigma_{23}')
view([0 90])
caxis([min_S 700])
axis equal tight

figure(7)
surf(InputData.y_grid, InputData.z_grid, SEFF)
title('\Sigma_{VM}')
view([0 90])
caxis([min_S 700])
axis equal tight