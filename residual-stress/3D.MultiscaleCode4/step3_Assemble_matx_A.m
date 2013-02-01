clear all
close all
clc

%%% INPUT DATA GENERATED IN STEP0
load inputdata.A1.mat

PFNAME_A        = fullfile(inputdata.PNAME_SOLUTION, inputdata.FNAME_A);
PFNAME_N        = fullfile(inputdata.PNAME_SOLUTION, inputdata.FNAME_N);
PFNAME_L        = fullfile(inputdata.PNAME_SOLUTION, inputdata.FNAME_L);
PFNAME_Nz       = fullfile(inputdata.PNAME_SOLUTION, inputdata.FNAME_Nz);
PFNAME_Nxy      = fullfile(inputdata.PNAME_SOLUTION, inputdata.FNAME_Nxy);
PFNAME_Fvec     = fullfile(inputdata.PNAME_SOLUTION, inputdata.FNAME_Fvec);
PFNAME_LS       = fullfile(inputdata.PNAME_SOLUTION, inputdata.FNAME_LS);
PFNAME_Fmk      = fullfile(inputdata.PNAME_SOLUTION, inputdata.FNAME_Fmk);
PFNAME_Re       = fullfile(inputdata.PNAME_SOLUTION, inputdata.FNAME_Re);
PFNAME_RANK     = fullfile(inputdata.PNAME_SOLUTION, inputdata.FNAME_RANK);
PFNAME_INVIG    = fullfile(inputdata.PNAME_SOLUTION, inputdata.FNAME_INVIG);
PFNAME_GRID     = fullfile(inputdata.PNAME_SOLUTION, inputdata.FNAME_GRID);

nx  = length(inputdata.x_grid);
ny  = length(inputdata.y_grid);
nz  = length(inputdata.z_grid);

Alpha   = [0 30 60 90];
na      = length(Alpha);

%%%
% load texture
pfname_odf  = fullfile(inputdata.pname_odf, inputdata.fname_odf);
odf_data    = load(pfname_odf);

%%%
% orientation space mesh
frmesh      = odf_data.wscub.frmesh;
nnode       = frmesh.numind;
sum_l2      = sum(frmesh.l2ip*odf_data.odf);
RMats       = RMatOfQuat(QuatOfRod(frmesh.crd(:,1:nnode)));

%%% GENERATE SH MODES
if inputdata.Use20xMesh
    load wscub20x
    [els, crds] = MeshCoordinates(wscub20x.frmesh, frmesh.crd(:,1:frmesh.numind));
    
    dh      = zeros(frmesh.numind,inputdata.nHarm);
    eVals   = zeros(inputdata.nHarm,1);
    for iHarm = 1:1:inputdata.nHarm
        dh(:,iHarm)     = EvalMeshFunc(wscub20x.frmesh, ...
            wscub20x.frmesh.dh(:,iHarm), els, crds);
        eVals(iHarm)    = wscub20x.frmesh.eVals(iHarm);
        
        %%% DH IN FRMESH
        [faces, multiplicity]   = MeshFaces(frmesh.con);
        scon                    = faces(:, (multiplicity == 1));
        surfmesh                = MeshStructure(frmesh.crd, scon, frmesh.eqv);
        
%         figure(1)
%         subplot(4,6,iHarm)
%         PlotSurface(surfmesh, dh(:,iHarm), ...
%             'ShowMesh',  'on');
        
%         PlotFR(frmesh, dh(:,iHarm), ...
%             'Clim', [-5 5]);
        
%         %%% ORIGINAL DH
%         [faces, multiplicity]   = MeshFaces(wscub20x.frmesh.con);
%         scon                    = faces(:, (multiplicity == 1));
%         surfmesh                = MeshStructure(wscub20x.frmesh.crd, scon, wscub20x.frmesh.eqv);
%         
%         figure(2)
%         subplot(4,6,iHarm)
%         PlotSurface(surfmesh, wscub20x.frmesh.dh(:,iHarm), ...
%             'ShowMesh',  'off');
    end
else
    [dh, eVals] = DiscreteHarmonics(frmesh, inputdata.nHarm);
end

%%%
% SX-STIFFNESS
C   = BuildElasticityMatrix(inputdata.C);
S   = BuildElasticityMatrix(inputdata.S);

%%%
% ASSEMBLE DH MATRIX
DH	= zeros(6*nnode,6*inputdata.nHarm);

DH(1:6:end-5,1:inputdata.nHarm)                     = dh;           % 11
DH(2:6:end-4,inputdata.nHarm+1:2*inputdata.nHarm)   = dh;           % 22
DH(3:6:end-3,2*inputdata.nHarm+1:3*inputdata.nHarm) = dh;           % 33
DH(4:6:end-2,3*inputdata.nHarm+1:4*inputdata.nHarm) = dh;           % 12
DH(5:6:end-1,4*inputdata.nHarm+1:5*inputdata.nHarm) = dh;           % 13
DH(6:6:end,5*inputdata.nHarm+1:6*inputdata.nHarm)   = dh;           % 23

DH	= sparse(DH);

%%%
% ASSEMBLE BigODF MATRIX
% NEEDED FOR SUMMATION OPERATOR N TO SUM STRESSES OVER OS
BigODF              = zeros(6*nnode,6);
BigODF(1:6:end-5,1) = odf_data.odf;
BigODF(2:6:end-4,2) = odf_data.odf;
BigODF(3:6:end-3,3) = odf_data.odf;
BigODF(4:6:end-2,4) = odf_data.odf;
BigODF(5:6:end-1,5) = odf_data.odf;
BigODF(6:6:end,6)   = odf_data.odf;
BigODF              = sparse(BigODF);

%%%
% BUILD BigC or BigS - STIFFNESS / COMPLIANCE MATRIX
if inputdata.useLSDF
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
% INITIALIZE VARIABLES
q   = [];
eqq = [];
eqq_mk  = [];

BigA        = [];
BigNMat     = [];
BigNMatz    = [];
BigNMatij   = [];
BigL        = [];
Fvec        = [];

PosNum  = 1;
q_last  = 0;
HarmNum     = 0;
SijNum      = 0;

% STRESSES
SIGMA11 = zeros(na,ny,nz);
SIGMA22 = zeros(na,ny,nz);
SIGMA33 = zeros(na,ny,nz);
SIGMA12 = zeros(na,ny,nz);
SIGMA13 = zeros(na,ny,nz);
SIGMA23 = zeros(na,ny,nz);
STRESS  = [];

% WEIGHTED ERROR
Re      = cell(size(inputdata.hkls,1)+1,1);
% RANK
Rk      = zeros(nx*ny*nz,4);

tic
ndv = 0;
IG  = [];
UB  = [];
LB  = [];
for i = 1:1:nx
    for j = 1:1:ny
        for ii = 1:1:na
            for k = 1:1:nz
                if ii == 1
                    input_fname = 'inputdata.A1.mat';
                elseif ii == 2
                    input_fname = 'inputdata.A2.mat';
                elseif ii == 3
                    input_fname = 'inputdata.A3.mat';
                elseif ii == 4
                    input_fname = 'inputdata.A4.mat';
                end
                inputdata_alpha = getfield(load(input_fname), 'inputdata');
                
                XSi = inputdata.x_grid(i);
                YSj = inputdata.y_grid(j);
                ZSk = inputdata.z_grid(k);
                a   = Alpha(ii);
                
                uvw = [...
                    'x=', num2str(XSi, '%10.3f'), '.', ...
                    'y=', num2str(YSj, '%10.3f'), '.', ...
                    'z=', num2str(ZSk, '%10.3f')
                    ];
                
                % LOAD SCATTERING VECTORS AND LATTICE STRAINS PER DV
                q_DV    = [];
                eqq_DV  = [];
                q_size_DV   = [];
                q_last_DV   = 0;
                for u = 1:1:size(inputdata_alpha.hkls,1)
                    H   = inputdata_alpha.hkls(u,1);
                    K   = inputdata_alpha.hkls(u,2);
                    L   = inputdata_alpha.hkls(u,3);
                    
                    f   = [num2str(H), num2str(K), num2str(L), '.LS.flt.srt.uniq.data'];
                    pf  = fullfile(inputdata_alpha.PNAME_LS_FILTER, uvw, f);
                    
                    disp(['pos num : ', num2str(PosNum), ' - loading ', f])
                    data    = load(pf);
                    
                    % UPDATE SCATTERING VECTORS AND LATTICE STRAINS PER DV
                    q_DV    = [q_DV; data(:,5:7)];
                    eqq_DV  = [eqq_DV; data(:,8)];
                    
                    idx_q_DV    = size(data,1) + q_last_DV;
                    q_size_DV   = [q_size_DV; ...
                        PosNum inputdata.hkls(u,:) idx_q_DV];
                    
                    q_last_DV   = q_size_DV(end,5);
                    
                    % FOR ALL DVs
                    q   = [q; data(:,1:3)];
                    eqq = [eqq; data(:,8)]; % RHS OF THE SOLVER
                    
                    idx_q   = size(data,1) + q_last;
                    eqq_mk  = [eqq_mk; ...
                        PosNum inputdata_alpha.hkls(u,:) idx_q];
                    
                    q_last  = eqq_mk(end,5);
                    
                    % figure(1)
                    % subplot(1,3,u)
                    % PlotSPF(data(:,5:7), data(:,8))
                end
                
                % TIME TO GENERATE INVERSION MATRIX
                % LOAD A MATRIX FOR DV : VARIABLE NAME IS A
                if inputdata.useLSDF
                    f   = ['A.DV.x=', num2str(XSi, '%10.3f'), ...
                        '.y=', num2str(YSj, '%10.3f'), ...
                        '.z=', num2str(ZSk, '%10.3f'), ...
                        '.mat'];
                else
                    f   = ['A.SODF2LS.DV.x=', num2str(XSi, '%10.3f'), ...
                        '.y=', num2str(YSj, '%10.3f'), ...
                        '.z=', num2str(ZSk, '%10.3f'), ...
                        '.mat'];
                end
                pf  = fullfile(inputdata_alpha.PNAME_LS_FILTER, f);
                load(pf)
                
                %%%
                % ASSEMBLE A INTO BigA
                A   = A*DH;
                if inputdata.UsePM
                    ATA     = A'*A;
                    
                    ri  = size(BigA,1) + 1;
                    rf  = size(BigA,1) + inputdata.nHarm*6;
                    ci  = inputdata.nHarm*6*(PosNum-1) + 1;
                    cf  = inputdata.nHarm*6*PosNum;
                    BigA(ri:rf, ci:cf)  = ATA;
                    
                    ri  = size(BigL,1) + 1;
                    rf  = size(BigL,1) + size(A,1);
                    ci  = inputdata.nHarm*6*(PosNum-1) + 1;
                    cf  = inputdata.nHarm*6*PosNum;
                    BigL(ri:rf, ci:cf)  = A;
                    
                    Fvec    = [Fvec; A'*eqq_DV];
                else
                    ri  = size(BigA,1) + 1;
                    rf  = size(BigA,1) + size(A,1);
                    ci  = inputdata.nHarm*6*(PosNum-1) + 1;
                    cf  = inputdata.nHarm*6*PosNum;
                    BigA(ri:rf, ci:cf)  = A;
                    
                    Fvec    = [Fvec; eqq_DV];
                end
                
                %%%
                % SUMMATION OPERATOR TO SUM STRESS OVER ORIENTATION
                if inputdata.useLSDF
                    Nmat	= BigODF'*BigL2IP*BigC*DH/sum_l2;
                else
                    Nmat	= BigODF'*BigL2IP*DH/sum_l2;
                end
                
                %%% GENERATE SIGMA_APPLIED TO BE USED FOR IG GENERATION
                %%% USE SIMPLE INVERSION TO GENERATE APPRAENT MACRO STRESS
                if inputdata.UsePM
                    X_DV        = ATA\(A'*eqq_DV);    %%% SIMPLEST INV
                else
                    X_DV        = A\eqq_DV;           %%% SIMPLEST INV
                end
                SIGMA_APPLIED       = Nmat*X_DV;
                SIGMA_APPLIED(4:6)  = [0 0 0]';
                
                %%%
                % PERFORM SIMPLE INVERSION
                %%%
                %%% FIRST GENERATE IG, UB, LB
                %%% CALCULATE INV-IG
                % SIGMA_APPLIED   = [1000 -1000 500 0 0 0]';
                S_MACRO         = AggregateElasticityMatrix(frmesh, odf_data.odf, S);
                EPSILON_APPLIED = S_MACRO*SIGMA_APPLIED;
                ig              = CalculateSPFInversionIG(frmesh, dh, ...
                    EPSILON_APPLIED, C);
                
                k2  = 0.5;
                ub  = ig + abs(ig)*k2;
                lb  = ig - abs(ig)*k2;
                
                iub = abs(ig) < 0.5;
                ilb = abs(ig) < 0.5;
                
                ub(iub) = +1.0;
                lb(ilb) = -1.0;
                
                IG  = [IG; ig];
                UB  = [UB; ub];
                LB  = [LB; lb];
                
                %%% PLOT IG
                SHIG    = reshape(ig', inputdata.nHarm, 6);
                SHIG_11 = SHIG(:,1);
                SHIG_22 = SHIG(:,2);
                SHIG_33 = SHIG(:,3);
                SHIG_12 = SHIG(:,4);
                SHIG_13 = SHIG(:,5);
                SHIG_23 = SHIG(:,6);
                
                SODF11_IG   = dh*SHIG_11;
                SODF22_IG   = dh*SHIG_22;
                SODF33_IG   = dh*SHIG_33;
                SODF12_IG   = dh*SHIG_12;
                SODF13_IG   = dh*SHIG_13;
                SODF23_IG   = dh*SHIG_23;
                
                % PlotFR(frmesh, SODF11_IG)
                % PlotFR(frmesh, SODF22_IG)
                % PlotFR(frmesh, SODF33_IG)
                % PlotFR(frmesh, SODF12_IG)
                % PlotFR(frmesh, SODF13_IG)
                % PlotFR(frmesh, SODF23_IG)
                
                %%% INVERT
                if inputdata.UsePM
                    % X_DV        = ATA\(A'*eqq_DV);    %%% SIMPLEST INV
                    X_DV    = lsqlin(ATA, A'*eqq_DV, ...
                        [], [], ...
                        [], [], ...
                        lb, ub, ...
                        ig);
                    eqq_calc_DV = A*X_DV;
                else
                    X_DV        = A\eqq_DV;           %%% SIMPLEST INV
                    X_DV    = lsqlin(A, eqq_DV, ...
                        [], [], ...
                        [], [], ...
                        lb, ub, ...
                        ig);
                    eqq_calc_DV = A*X_DV;
                end
                
                %%% PLOT INV RESULT
                SODF    = DH*X_DV;
                SODF11  = SODF(1:6:end-5);
                SODF22  = SODF(2:6:end-4);
                SODF33  = SODF(3:6:end-3);
                SODF12  = SODF(4:6:end-2);
                SODF13  = SODF(5:6:end-1);
                SODF23  = SODF(6:6:end-0);
                
%                 PlotFR(frmesh, SODF11)
%                 PlotFR(frmesh, SODF22)
%                 PlotFR(frmesh, SODF33)
%                 PlotFR(frmesh, SODF12)
%                 PlotFR(frmesh, SODF13)
%                 PlotFR(frmesh, SODF23)
                
                SIGMA	= Nmat*X_DV;        % MPa (DEPENDS ON C INPUT)
                disp('radpos alpha Sxx Syy Szz Syz Sxz Sxy FOR IG')
                disp([j a round(SIGMA_APPLIED')])
                disp('radpos alpha Sxx Syy Szz Syz Sxz Sxy SOL')
                disp([j a round(SIGMA')])
                disp('% mean LS diff')
                MEAN_LS_DIFF(PosNum,1)  = mean(abs(eqq_calc_DV - eqq_DV))/max(abs(eqq_DV));
                disp(mean(abs(eqq_calc_DV - eqq_DV)))
                disp('% max LS diff')
                MAX_LS_DIFF(PosNum,1)   = max(abs(eqq_calc_DV - eqq_DV))/max(abs(eqq_DV));
                disp(max(abs(eqq_calc_DV - eqq_DV)))
                
                SIGMA11(ii,j,k)     = SIGMA(1,1);
                SIGMA22(ii,j,k)     = SIGMA(2,1);
                SIGMA33(ii,j,k)     = SIGMA(3,1);
                SIGMA12(ii,j,k)     = SIGMA(6,1);
                SIGMA13(ii,j,k)     = SIGMA(5,1);
                SIGMA23(ii,j,k)     = SIGMA(4,1);
                STRESS(PosNum,:)    = SIGMA';
                
                % TO PASS CRYSTAL AVE STRESS TO MACRO CODE
                % SIGMA NEEDS TO TO BE TRANSFORMED TO SAMPLE COORDINATE SYSTEM
                R   = [
                    cosd(a) -sind(a) 0;
                    sind(a)  cosd(a) 0;
                    0 0 1;
                    ];
                
                xyz0    = [XSi -YSj+6.35 ZSk]';
                xyz     = R*xyz0;
                T       = VectorizedCOBMatrix(R);
                Nmat    = T*Nmat;
                
                % figure(1000)
                % axis([-7 7 -7 7 -1 1])
                % view([-35 68])
                % hold on
                % grid on
                % plot3(xyz(1), xyz(2), xyz(3), 'k.')
                
                ndv = ndv + 1;
                dvc_x(1,ndv)    = xyz(1);
                dvc_y(1,ndv)    = xyz(2);
                dvc_z(1,ndv)    = xyz(3);
                
                %%% 3D MACRO CODE TAKES SHEARS W/O sqrt(2)
                %%% 3D MACRO CODE STRESS ORDER 11 22 33 23 13 12
                Nmat([4 5 6],:) = Nmat([6 5 4],:)./sqrt(2);
                
                ri  = SijNum+1;
                rf  = SijNum+6;
                ci  = HarmNum+1;
                cf  = HarmNum+6*inputdata.nHarm;
                BigNMat(ri:rf,ci:cf)    = Nmat;
                
                ri  = SijNum/2+1;
                rf  = SijNum/2+3;
                ci  = HarmNum+1;
                cf  = HarmNum+6*inputdata.nHarm;
                BigNMatz(ri:rf,ci:cf)   = Nmat([3 5 6],:);
                BigNMatij(ri:rf,ci:cf)  = Nmat([1 2 4],:);
                
                %%%
                % MARGIN OF ERROR CALCULATION
                % SPF IMAGE GENERATION
                Re{end,1}   = [Re{end,1}; ...
                    inputdata.x_grid(i), ...
                    inputdata.y_grid(j), ...
                    inputdata.z_grid(k), ...
                    ErrorRe(eqq_DV, eqq_calc_DV, 'Threshold', 0.0001)];
                Rk(PosNum,:)    = [inputdata.x_grid(i), ...
                    inputdata.y_grid(j), ...
                    inputdata.z_grid(k), ...
                    rank(full(A))/size(A,2)];
                for u = 1:1:size(inputdata.hkls,1)
                    H   = inputdata_alpha.hkls(u,1);
                    K   = inputdata_alpha.hkls(u,2);
                    L   = inputdata_alpha.hkls(u,3);
                    
                    if u == 1
                        ri  = 1;
                        rf  = q_size_DV(u,5);
                    else
                        ri  = q_size_DV(u-1,5) + 1;
                        rf  = q_size_DV(u,5);
                    end
                    recalc_data = [q_DV(ri:rf,:) eqq_calc_DV(ri:rf)];
                    org_data    = [q_DV(ri:rf,:) eqq_DV(ri:rf)];
                    
                    % CALCULATE ERROR MEASURES
                    Re{u,1} = [Re{u,1}; ...
                        inputdata.x_grid(i), ...
                        inputdata.y_grid(j), ...
                        inputdata.z_grid(k), ...
                        ErrorRe(eqq_DV(ri:rf), eqq_calc_DV(ri:rf),'Threshold', 0.0001)];
                    
                    % SAVE RECALCULATED DATA
                    if inputdata.eqq_calc_fig
                        if inputdata.useLSDF
                            f   = [num2str(H), num2str(K), num2str(L), '.LS.recalc.data'];
                        else
                            f   = [num2str(H), num2str(K), num2str(L), '.LS.recalc.SODF2LS.data'];
                        end
                        p   = fullfile(inputdata_alpha.PNAME_LS_FILTER, uvw);
                        pf  = fullfile(p,f);
                        
                        % save recalculated data
                        disp(['pos num : ', num2str(PosNum), ' - saving ', f])
                        fid = fopen(pf, 'w');
                        fprintf(fid, '%% RECALCULATED SPF DATA VIA INVERSION W/O MACRO BCS\n');
                        fprintf(fid, '%% q eqq\n');
                        fprintf(fid, '%2.6f %2.6f %2.6f %2.6f \n', recalc_data');
                        fclose(fid);
                        
                        % save recalculated data spf
%                         fig1    = figure(1000);
%                         PlotSPF(recalc_data(:,1:3), recalc_data(:,4), ...
%                             'Title', 'SPF_{CALC}', ...
%                             'DataRange', [-0.0056 0.0074])
%                         hold off
                        
                        if inputdata.useLSDF
                            f   = [num2str(H), num2str(K), num2str(L), '.SPF.recalc'];
                        else
                            f   = [num2str(H), num2str(K), num2str(L), '.SPF.recalc.SODF2LS'];
                        end
                        p   = fullfile(inputdata_alpha.PNAME_LS_FILTER, uvw);
                        pf  = fullfile(p,f);
                        
                        disp(['pos num : ', num2str(PosNum), ' - saving ', f, '.png'])
                        %saveas(fig1, [pf, '.png'], 'png')
                        disp(['pos num : ', num2str(PosNum), ' - saving ', f, '.fig'])
                        %saveas(fig1, [pf, '.fig'], 'fig')
                        
                        % save data - recalculated data spf
%                         fig2    = figure(2000);
%                         PlotSPF(recalc_data(:,1:3), org_data(:,4) - recalc_data(:,4), ...
%                             'Title', 'SPF_{EXP} - SPF_{CALC}')
%                         hold off
                        
                        if inputdata.useLSDF
                            f   = [num2str(H), num2str(K), num2str(L), '.delta_SPF'];
                        else
                            f   = [num2str(H), num2str(K), num2str(L), '.delta_SPF.SODF2LS'];
                        end
                        p   = fullfile(inputdata_alpha.PNAME_LS_FILTER, uvw);
                        pf  = fullfile(p,f);
                        
                        % save data - difference spf
                        disp(['pos num : ', num2str(PosNum), ' - saving ', f, '.png'])
%                         saveas(fig2, [pf, '.png'], 'png')
                        disp(['pos num : ', num2str(PosNum), ' - saving ', f, '.fig'])
%                         saveas(fig2, [pf, '.fig'], 'fig')
                        
%                         fig3    = figure(3000);
%                         PlotSPF(recalc_data(:,1:3), org_data(:,4), ...
%                             'Title', 'SPF_{EXP}', ...
%                             'DataRange', [-0.0056 0.0074])
%                         hold off
                    else
                        disp('not saving recalculated strain data ...')
                    end
                    pause(1)
                end
                % pause
                close all
                
                % UPDATE COUNTERS
                HarmNum = HarmNum + 6*inputdata.nHarm;
                SijNum  = SijNum + 6;
                PosNum  = PosNum + 1;
                
                % MEMORY MANAGEMENT
                clear A A2 Nmat
                BigA        = sparse(BigA);
                BigNMat     = sparse(BigNMat);
                BigNMatz    = sparse(BigNMatz);
                BigNMatij   = sparse(BigNMatij);
                if inputdata.UsePM
                    BigL    = sparse(BigL);
                end
            end
        end
    end
end

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

disp('generating DV grid file')
save(PFNAME_GRID, 'dvc_x', 'dvc_y', 'dvc_z', 'ndv')

%%%
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
surf(inputdata.y_grid, inputdata.z_grid, squeeze(S11))
title('\Sigma_{11}')
view([0 90])
caxis([min_S 700])
axis equal tight

figure(2)
surf(inputdata.y_grid, inputdata.z_grid, S22)
title('\Sigma_{22}')
view([0 90])
caxis([min_S 700])
axis equal tight

figure(3)
surf(inputdata.y_grid, inputdata.z_grid, S33)
title('\Sigma_{33}')
view([0 90])
caxis([min_S 700])
axis equal tight

figure(4)
surf(inputdata.y_grid, inputdata.z_grid, S12)
title('\Sigma_{12}')
view([0 90])
caxis([min_S 700])
axis equal tight

figure(5)
surf(inputdata.y_grid, inputdata.z_grid, S13)
title('\Sigma_{13}')
view([0 90])
caxis([min_S 700])
axis equal tight

figure(6)
surf(inputdata.y_grid, inputdata.z_grid, S23)
title('\Sigma_{23}')
view([0 90])
caxis([min_S 700])
axis equal tight

figure(7)
surf(inputdata.y_grid, inputdata.z_grid, SEFF)
title('\Sigma_{VM}')
view([0 90])
caxis([min_S 700])
axis equal tight