clear all
close all
clc

load inputdata.A1.mat
load Residual_Stress_3D_MLSEL_v1/meshdata.mat
load('/home/jspark/Documents/work/prjResidualStress/3D.MultiscaleCode3/A1234.SH.MLS.0,0,0,6.mat')
load('/home/jspark/Documents/work/Data/APS201103/GE/A1234.INVIG.mat')

Alpha       = [0 30 60 90];

x_target    = 0;
y_target    = -0.900;
z_target    = 0.2250;
a_target    = 90;

i_target    = find(inputdata.x_grid == x_target);
j_target    = find(inputdata.y_grid == y_target);
k_target    = find(inputdata.z_grid == z_target);
ii_target   = find(Alpha == a_target);

PFNAME_A    = fullfile(inputdata.PNAME_SAM, inputdata.FNAME_A);
PFNAME_N    = fullfile(inputdata.PNAME_SAM, inputdata.FNAME_N);
PFNAME_L    = fullfile(inputdata.PNAME_SAM, inputdata.FNAME_L);
PFNAME_Nz   = fullfile(inputdata.PNAME_SAM, inputdata.FNAME_Nz);
PFNAME_Nxy  = fullfile(inputdata.PNAME_SAM, inputdata.FNAME_Nxy);
PFNAME_Fvec = fullfile(inputdata.PNAME_SAM, inputdata.FNAME_Fvec);
PFNAME_LS   = fullfile(inputdata.PNAME_SAM, inputdata.FNAME_LS);
PFNAME_Re   = fullfile(inputdata.PNAME_SAM, inputdata.FNAME_Re);
PFNAME_RANK = fullfile(inputdata.PNAME_SAM, inputdata.FNAME_RANK);

load(PFNAME_LS)

nx  = length(inputdata.x_grid);
ny  = length(inputdata.y_grid);
nz  = length(inputdata.z_grid);
na  = length(Alpha);

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
    end
else
    [dh, eVals] = DiscreteHarmonics(frmesh, inputdata.nHarm);
end

%%%
% ASSEMBLE DH MATRIX
DH	= zeros(6*nnode,6*inputdata.nHarm);
DH(1:6:end-5,1:inputdata.nHarm)                     = dh;           % 11
DH(2:6:end-4,inputdata.nHarm+1:2*inputdata.nHarm)   = dh;           % 22
DH(3:6:end-3,2*inputdata.nHarm+1:3*inputdata.nHarm) = dh;           % 33
DH(4:6:end-2,3*inputdata.nHarm+1:4*inputdata.nHarm) = dh;   % 12
DH(5:6:end-1,4*inputdata.nHarm+1:5*inputdata.nHarm) = dh;   % 13
DH(6:6:end,5*inputdata.nHarm+1:6*inputdata.nHarm)   = dh;   % 23
DH	= sparse(DH);

npf = size(inputdata.hkls,1);
X   = SH(numnp*6+1:end);

load wscub3x
load wscub20x
load SDF_LSHR-Triax3.mat

dh1 = dh;
dh2 = zeros(wscub3x.frmesh.numind, inputdata.nHarm);
for iHarm = 1:1:inputdata.nHarm
    [els, crds]     = MeshCoordinates(wscub20x.frmesh, ...
        wscub3x.frmesh.crd(:,1:wscub3x.frmesh.numind));
    
    dh2(:,iHarm)    = EvalMeshFunc(wscub20x.frmesh, ...
        wscub20x.frmesh.dh(:,iHarm), els, crds);
    
    %%% DH IN FRMESH
    [faces, multiplicity]   = MeshFaces(wscub3x.frmesh.con);
    scon                    = faces(:, (multiplicity == 1));
    surfmesh                = MeshStructure(wscub3x.frmesh.crd, scon, wscub3x.frmesh.eqv);
    
%     figure(2)
%     subplot(4,6,iHarm)
%     PlotSurface(surfmesh, dh2(:,iHarm), ...
%         'ShowMesh',  'on');
end

c2(:,1) = GetSHCoeffs(sdf(1,:), wscub3x.frmesh, dh2);
c2(:,2) = GetSHCoeffs(sdf(4,:), wscub3x.frmesh, dh2);
c2(:,3) = GetSHCoeffs(sdf(6,:), wscub3x.frmesh, dh2);
c2(:,4) = GetSHCoeffs(sdf(2,:), wscub3x.frmesh, dh2);
c2(:,5) = GetSHCoeffs(sdf(3,:), wscub3x.frmesh, dh2);
c2(:,6) = GetSHCoeffs(sdf(5,:), wscub3x.frmesh, dh2);

SHC11	= zeros(nx*ny*na*nz,inputdata.nHarm);
SHC22	= zeros(nx*ny*na*nz,inputdata.nHarm);
SHC33	= zeros(nx*ny*na*nz,inputdata.nHarm);
SHC12	= zeros(nx*ny*na*nz,inputdata.nHarm);
SHC13	= zeros(nx*ny*na*nz,inputdata.nHarm);
SHC23   = zeros(nx*ny*na*nz,inputdata.nHarm);
ct  = 1;
delta_LS    = [];
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
                
                % LOAD A MATRIX FOR DV : VARIABLE NAME IS A
                f   = ['A.DV.x=', num2str(XSi, '%10.3f'), ...
                    '.y=', num2str(YSj, '%10.3f'), ...
                    '.z=', num2str(ZSk, '%10.3f'), ...
                    '.mat'];
                pf  = fullfile(inputdata_alpha.PNAME_LS_FILTER, f);
                disp(pf)
                load(pf)
                
                ri  = (ct-1)*6*inputdata.nHarm + 1;
                rf  = ct*6*inputdata.nHarm;
                X_DV    = X(ri:rf,1);
                IG_DV   = IG(ri:rf,1);
                UB_DV   = UB(ri:rf,1);
                LB_DV   = LB(ri:rf,1);
                
                SODF    = DH*X_DV;
                SODF_IG = DH*IG_DV;
                SODF_UB = DH*UB_DV;
                SODF_LB = DH*LB_DV;
                
                %%% EVALULATE SH COEFFICIENTS
                c1(:,1)     = GetSHCoeffs(SODF(1:6:end-5)', frmesh, dh1);
                c1(:,2)     = GetSHCoeffs(SODF(2:6:end-4)', frmesh, dh1);
                c1(:,3)     = GetSHCoeffs(SODF(3:6:end-3)', frmesh, dh1);
                c1(:,4)     = GetSHCoeffs(SODF(4:6:end-2)', frmesh, dh1)./sqrt(2);
                c1(:,5)     = GetSHCoeffs(SODF(5:6:end-1)', frmesh, dh1)./sqrt(2);
                c1(:,6)     = GetSHCoeffs(SODF(6:6:end-0)', frmesh, dh1)./sqrt(2);
                
                cIG(:,1)    = GetSHCoeffs(SODF_IG(1:6:end-5)', frmesh, dh1);
                cIG(:,2)    = GetSHCoeffs(SODF_IG(2:6:end-4)', frmesh, dh1);
                cIG(:,3)    = GetSHCoeffs(SODF_IG(3:6:end-3)', frmesh, dh1);
                cIG(:,4)    = GetSHCoeffs(SODF_IG(4:6:end-2)', frmesh, dh1)./sqrt(2);
                cIG(:,5)    = GetSHCoeffs(SODF_IG(5:6:end-1)', frmesh, dh1)./sqrt(2);
                cIG(:,6)    = GetSHCoeffs(SODF_IG(6:6:end-0)', frmesh, dh1)./sqrt(2);
                
                cUB(:,1)    = GetSHCoeffs(SODF_UB(1:6:end-5)', frmesh, dh1);
                cUB(:,2)    = GetSHCoeffs(SODF_UB(2:6:end-4)', frmesh, dh1);
                cUB(:,3)    = GetSHCoeffs(SODF_UB(3:6:end-3)', frmesh, dh1);
                cUB(:,4)    = GetSHCoeffs(SODF_UB(4:6:end-2)', frmesh, dh1)./sqrt(2);
                cUB(:,5)    = GetSHCoeffs(SODF_UB(5:6:end-1)', frmesh, dh1)./sqrt(2);
                cUB(:,6)    = GetSHCoeffs(SODF_UB(6:6:end-0)', frmesh, dh1)./sqrt(2);
                
                cLB(:,1)    = GetSHCoeffs(SODF_LB(1:6:end-5)', frmesh, dh1);
                cLB(:,2)    = GetSHCoeffs(SODF_LB(2:6:end-4)', frmesh, dh1);
                cLB(:,3)    = GetSHCoeffs(SODF_LB(3:6:end-3)', frmesh, dh1);
                cLB(:,4)    = GetSHCoeffs(SODF_LB(4:6:end-2)', frmesh, dh1)./sqrt(2);
                cLB(:,5)    = GetSHCoeffs(SODF_LB(5:6:end-1)', frmesh, dh1)./sqrt(2);
                cLB(:,6)    = GetSHCoeffs(SODF_LB(6:6:end-0)', frmesh, dh1)./sqrt(2);
                
                SHC11(ct,:) = c1(:,1);
                SHC22(ct,:) = c1(:,2);
                SHC33(ct,:) = c1(:,3);
                SHC12(ct,:) = c1(:,4);
                SHC13(ct,:) = c1(:,5);
                SHC23(ct,:) = c1(:,6);
                
                %%% PLOTTING FOR A PARTICULAR DV STARTS HERE
                if i == i_target && ...
                        j == j_target && ...
                        ii == ii_target && ...
                        k == k_target
                    disp('plotting simulation sodf as a reference')
                    PlotFR(wscub3x.frmesh, sdf(1,:), 'CameraPosition', [3.7, 5.5, 2.6])
                    PlotFR(wscub3x.frmesh, sdf(4,:), 'CameraPosition', [3.7, 5.5, 2.6])
                    PlotFR(wscub3x.frmesh, sdf(6,:), 'CameraPosition', [3.7, 5.5, 2.6])
                    PlotFR(wscub3x.frmesh, sdf(2,:), 'CameraPosition', [3.7, 5.5, 2.6])
                    PlotFR(wscub3x.frmesh, sdf(3,:), 'CameraPosition', [3.7, 5.5, 2.6])
                    PlotFR(wscub3x.frmesh, sdf(5,:), 'CameraPosition', [3.7, 5.5, 2.6])
                    
                    %%%
                    disp('plotting experimental sodf at a particular DV ...')
                    PlotFR(frmesh, SODF(1:6:end-5), 'CameraPosition', [3.7, 5.5, 2.6])
                    PlotFR(frmesh, SODF(2:6:end-4), 'CameraPosition', [3.7, 5.5, 2.6])
                    PlotFR(frmesh, SODF(3:6:end-3), 'CameraPosition', [3.7, 5.5, 2.6])
                    PlotFR(frmesh, SODF(4:6:end-2)./sqrt(2), 'CameraPosition', [3.7, 5.5, 2.6])
                    PlotFR(frmesh, SODF(5:6:end-1)./sqrt(2), 'CameraPosition', [3.7, 5.5, 2.6])
                    PlotFR(frmesh, SODF(6:6:end-0)./sqrt(2), 'CameraPosition', [3.7, 5.5, 2.6])
                    
                    disp(round([ ...
                        (SODF(1:6:end-5)'*frmesh.l2ip*odf_data.odf)./sum(frmesh.l2ip*odf_data.odf); ...
                        (SODF(2:6:end-4)'*frmesh.l2ip*odf_data.odf)./sum(frmesh.l2ip*odf_data.odf); ...
                        (SODF(3:6:end-3)'*frmesh.l2ip*odf_data.odf)./sum(frmesh.l2ip*odf_data.odf); ...
                        (SODF(4:6:end-2)'*frmesh.l2ip*odf_data.odf)./sum(frmesh.l2ip*odf_data.odf)./sqrt(2); ...
                        (SODF(5:6:end-1)'*frmesh.l2ip*odf_data.odf)./sum(frmesh.l2ip*odf_data.odf)./sqrt(2); ...
                        (SODF(6:6:end-0)'*frmesh.l2ip*odf_data.odf)./sum(frmesh.l2ip*odf_data.odf)./sqrt(2); ...
                        ]'))
                    
                    figure(5000)
                    xHarm   = 1:1:inputdata.nHarm;
                    for iComp = 1:1:6
                        subplot(2,3,iComp)
                        plot(xHarm, c1(1:inputdata.nHarm,iComp), 'r.-')
                        hold on
                        plot(xHarm, c2(1:inputdata.nHarm,iComp), 'b.-')
                        % plot(xHarm, cIG(1:inputdata.nHarm,iComp), 'g.-')
                        % plot(xHarm, cUB(1:inputdata.nHarm,iComp), 'm^-')
                        % plot(xHarm, cLB(1:inputdata.nHarm,iComp), 'mv-')
                        xlabel('SH Mode Number')
                        ylabel('SH Coefficient (MPa)')
                    end
                    hold off
                    return
                end
                close all
                ct  = ct + 1;
            end
        end
    end
end

figure(1)
subplot(2,3,1)
imagesc(SHC11)
colorbar vert
axis square
title('xx')

subplot(2,3,2)
imagesc(SHC22)
colorbar vert
axis square
title('yy')

subplot(2,3,3)
imagesc(SHC33)
colorbar vert
axis square
title('zz')

subplot(2,3,4)
imagesc(SHC12)
colorbar vert
axis square
title('xy')

subplot(2,3,5)
imagesc(SHC13)
colorbar vert
axis square
title('xz')

subplot(2,3,6)
imagesc(SHC23)
colorbar vert
axis square
title('yz')