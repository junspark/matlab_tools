clear all
close all
clc

tic
%%% INPUT DATA GENERATED IN STEP0
load inputdata.A1.mat

%%%
% load texture
pfname_odf  = fullfile(inputdata.pname_odf, inputdata.fname_odf);
odf_data    = load(pfname_odf);

%%%
% orientation space mesh
frmesh      = odf_data.wscub.frmesh;
nnode       = frmesh.numind;

%%%
% fiber is a tube in orientation space
% not really necessary for now ...
% we will have to look into this
pts=[1 0 0]';
[TH,PHI,R]  = cart2sph(pts(1,1),pts(2,1),pts(3,1));
V1  = deg2rad(2.5);
V2  = deg2rad(2.5)*cosd(45);
THpts   = [TH+V1; TH+V2; TH; TH-V2; TH-V1; TH-V2; TH; TH+V2];
PHIpts  = [PHI; PHI+V2; PHI+V1; PHI+V2; PHI; PHI-V2; PHI-V1; PHI-V2];
Rpts    = ones(size(THpts));

[x, y, z]   = sph2cart(THpts,PHIpts,Rpts);
pts         = [pts [x y z]'];
pts         = UnitVector(pts);

sph.crd = pts;
sph.con = [1 2 3; 1 3 4; 1 4 5; 1 5 6; 1 6 7; 1 7 8; 1 8 9; 1 9 2]';
sph.eqv = [];
qrule2d = LoadQuadrature('qr_trid06p12');
[gqr, sph.l2ip] = SphGQRule(sph, qrule2d);
w = sum(sph.l2ip)./sum(sum(sph.l2ip));

%%% BRAGG ANGLE CALCULATION
[~, th] = PlaneSpacings(inputdata.lattparms, ...
    'cubic', inputdata.hkls', ...
    inputdata.wavelength);

q       = [];
hkl_qq  = [];
for i = 1:1:size(inputdata.hkls,1)
    PHI = -deg2rad(th(i));
    PHI	= repmat(PHI, inputdata.eta_bins, 1);    
    R	= ones(inputdata.eta_bins, 1);
    
    THETA   = deg2rad(inputdata.eta_grid');
    [qbase(:,1), qbase(:,2) qbase(:,3)] = sph2cart(THETA, PHI, R);
    
    for j = 1:1:inputdata.omega_bins
        %%% SCATTERING VECTOR CALCULATION
        omega   = inputdata.omega_grid(j);
        Ry  = [...
            cosd(omega) 0 sind(omega); ...
            0 1 0; ...
            -sind(omega) 0 cosd(omega); ...
            ];
        qj  = Ry*qbase';
        q   = [q; qj'];
        
        hkl_qq  = [hkl_qq; ...
            ones(inputdata.eta_bins,1)*inputdata.hkls(i,:)];
    end
end

if inputdata.useLSDF
    disp('Generate LSDF <> LS ...')
    f   = 'BaseA.DV.mat';
    pf  = fullfile(inputdata.PNAME_LS_FILTER, f);
else
    disp('Generate LSDF <> LS ...')
    f   = 'BaseA.SODF2LS.DV.mat';
    pf  = fullfile(inputdata.PNAME_LS_FILTER, f);
end
q       = q';
hkl_qq  = hkl_qq';

A0  = sparse(size(q,2),6*nnode);
disp('generating Base A matrix ...')
tic
for ii = 1:1:size(q,2)
    if ~rem(ii,500)
        disp([num2str(ii/size(q,2)*100) '%'])
    end
    
    % STRAIN = {e11 e22 e33 sqrt(2)e12 sqrt(2)e13 sqrt(2)e23}
    % STRESS = {s11 s22 s33 sqrt(2)s12 sqrt(2)s13 sqrt(2)s23}
    if inputdata.useLSDF
        A0(ii,:)    = BuildLsdfSpfMatrix(hkl_qq(:,ii), ...
            frmesh, inputdata.sym, q(:,ii), inputdata.div, odf_data.odf, 1);    %%% FOR NOW TUBE WEIGHTING IS DISABLED
    else
        S   = BuildElasticityMatrix(inputdata.S);
        A0(ii,:)    = BuildSodfSpfMatrix(hkl_qq(:,ii), ...
            frmesh, inputdata.sym, q(:,ii), inputdata.div, odf_data.odf, 1, S); %%% FOR NOW TUBE WEIGHTING IS DISABLED
    end
end

toc
disp(['saving ', f])
save(pf, 'A0');

data    = [hkl_qq; q];
f   = 'BaseA.hkl.q.DV.data';
pf  = fullfile(inputdata.PNAME_LS_FILTER, f);
fid = fopen(pf, 'w');
fprintf(fid, '%% h k l qx qy qz\n');
fprintf(fid, '%f %f %f %f %f %f \n', data);
fclose(fid);