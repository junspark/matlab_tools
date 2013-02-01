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

p   = inputdata.PNAME_LS_FILTER;
%%%
% BUILD A MATRIX FOR EACH DV
for i = 1:1:length(inputdata.x_grid)
    for j = 1:1:length(inputdata.y_grid)
        for k = 1:1:length(inputdata.z_grid)
            tic
            q   = [];
            eqq = [];
            hkl_qq  = [];
            
            XSi = inputdata.x_grid(i);
            YSj = inputdata.y_grid(j);
            ZSk = inputdata.z_grid(k);
            
            uvw = [...
                'x=', num2str(XSi, '%10.3f'), '.', ...
                'y=', num2str(YSj, '%10.3f'), '.', ...
                'z=', num2str(ZSk, '%10.3f')
                ];
            
            for u = 1:1:size(inputdata.hkls,1)
                H   = inputdata.hkls(u,1);
                K   = inputdata.hkls(u,2);
                L   = inputdata.hkls(u,3);
                
                f   = [num2str(H), num2str(K), num2str(L), '.LS.flt.srt.uniq.data'];
                pf  = fullfile(p, uvw, f);
                
                disp(['loading ', f])
                data    = load(pf);
                
                q   = [q; data(:,5:7)];
                eqq = [eqq; data(:,8)];
                if size(data,1) > 1800
                    disp(['lattice strains more than 1800 ', num2str(size(data,1))])
                    return
                end
                hkl_qq  = [hkl_qq; ...
                    ones(size(data,1),1)*inputdata.hkls(u,:)];
            end
            % DV DEPENDENT TEXTURE HERE IF NEEDED
            % FOR NOW FOR EFFICIENCY TEXTURE IS LOADED OUTSIDE THE LOOP
            % AND ASSUMED UNIFORM TEXTURE
            
            % TIME TO GENERATE INVERSION MATRIX
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
            pf  = fullfile(inputdata.PNAME_LS_FILTER, f);
            
            q       = q';
            hkl_qq  = hkl_qq';
            
            A	= zeros(size(q,2),6*nnode);
            disp(['generating A matrix for ', f])
            if inputdata.useA0
                if inputdata.useLSDF
                    disp('A is LSDF <> LS ...')
                    f_BaseA     = 'BaseA.DV.mat';
                    pf_BaseA    = fullfile(inputdata.PNAME_LS_FILTER, f_BaseA);
                else
                    disp('A is SODF <> LS ...')
                    f_BaseA     = 'BaseA.SODF2LS.DV.mat';
                    pf_BaseA    = fullfile(inputdata.PNAME_LS_FILTER, f_BaseA);
                end
                load(pf_BaseA)
                
                f_BaseA     = 'BaseA.hkl.q.DV.data';
                pf_BaseA    = fullfile(inputdata.PNAME_LS_FILTER, f_BaseA);
                data        = load(pf_BaseA);
                hkl_qq_A0   = data(:,1:3)';
                q_A0        = data(:,4:6)';
                for ii = 1:size(q,2)
                    if ~rem(ii,500)
                        disp([num2str(ii/size(q,2)*100) '%'])
                    end
                    
                    ih  = hkl_qq(1,ii) == hkl_qq_A0(1,:);
                    ik  = hkl_qq(2,ii) == hkl_qq_A0(2,:);
                    il  = hkl_qq(3,ii) == hkl_qq_A0(3,:);
                    
                    idx = ih & ik & il;
                    
                    dqx = abs(q(1,ii) - q_A0(1,idx));
                    dqy = abs(q(2,ii) - q_A0(2,idx));
                    dqz = abs(q(3,ii) - q_A0(3,idx));
                    
                    iqx = dqx <= 0.0001;
                    iqy	= dqy <= 0.0001;
                    iqz	= dqz <= 0.0001;
                    
                    idx = iqx & iqy & iqz;
                    
                    if ~sum(idx)
                        disp('oh no! matching q not found ...')
                        return
                    end
                    if sum(idx) > 1
                        disp('oh no! multiple q match...')
                        return
                    end
                    A(ii,:) = A0(idx,:);
                end
            else
                for ii = 1:size(q,2)
                    if ~rem(ii,500)
                        disp([num2str(ii/size(q,2)*100) '%'])
                    end
                    
                    % STRAIN = {e11 e22 e33 sqrt(2)e12 sqrt(2)e13 sqrt(2)e23}
                    % STRESS = {s11 s22 s33 sqrt(2)s12 sqrt(2)s13 sqrt(2)s23}
                    if inputdata.useLSDF
                        A(ii,:) = BuildLsdfSpfMatrix(hkl_qq(:,ii), ...
                            frmesh, inputdata.sym, q(:,ii), inputdata.div, odf_data.odf, w);
                    else
                        A(ii,:) = BuildSodfSpfMatrix(hkl_qq(:,ii), ...
                            frmesh, inputdata.sym, q(:,ii), inputdata.div, odf_data.odf, w, S);
                    end
                end
            end
            A   = sparse(A);
            max(max(A(:)))
            if inputdata.saveA
                disp(['saving ', f])
                save(pf, 'A');
            else
                disp(['not saving ', f])
            end
            clear A
            toc
        end
    end
end