clear all
close all
clc

tic
%%% INPUT DATA GENERATED IN STEP0
InputData   = load('example_rsinput.mat');

%%%
% FIBER IS A TUBE IN ORIENTATION SPACE DUE TO OMEGA AND ETA STEP SIZES
% COMPUTE WEIGHTS FOR FIBER TUBE
w = FiberTubeWeight(InputData.deta, InputData.domega);

%%%
% PREP FRMESH
frmesh  = InputData.frmesh;
nnode   = frmesh.numind;

%%% GENERATE A MATRIX FOR EACH DV
ct  = 0;
for i = 1:1:(InputData.ndv/4)   % WE DUPLICATED DVS 
    %%% LOAD ODF
    disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
    disp('loading odf for dv ...')
    odf_data    = load(InputData.PFNAME_ODF{i,1});
    % PlotFR(frmesh, odf_data, 'Symmetries', 'hexagonal')
    
    %%% PREPARE SCATTERING VECTOR / HKL LIST
    disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
    q_DV    = [];
    hkl_DV  = [];
    for j = 1:1:InputData.nhkls
        fname_LS_DATA   = [ ...
            InputData.FNAME_LS_DATA_ROOT, ...
            num2str(InputData.hkls(j,1)), ...
            num2str(InputData.hkls(j,2)), ...
            num2str(InputData.hkls(j,3)), ...
            '.data'];
        disp('loading ...')
        disp(fname_LS_DATA)
        
        pfname_LS_DATA  = fullfile(InputData.PNAME_LS_DATA{i,1}, fname_LS_DATA);
        
        LS_data	= load(pfname_LS_DATA);
        
        %%% NEED TO CONVERT TO ORTHOGONAL BASIS IF HCP
        aspect  = InputData.LatticePrms(2)/InputData.LatticePrms(1);
        hkl     = ConvertMillerBravais(InputData.hkls(j,:)');
        hkl     = MillerBravaisToNormal(hkl, aspect)';
        hkl     = ones(size(LS_data,1),1)*hkl;
        
        hkl_DV  = [hkl_DV; hkl];
        q_DV    = [q_DV; LS_data(:,1:3)];
    end
    q_DV    = q_DV';
    hkl_DV  = hkl_DV';
    
    %%% START A MATRIX GENERATION
    disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
    disp('generating Base A matrix ...')
    tic
    
    if InputData.UseLSDF
        disp('Generate LSDF <> LS ...')
        f   = ['BaseA.LSDF2LS.DV', num2str(ct), '.mat'];
        pf  = fullfile(InputData.PNAME_SOLUTION, f);
    else
        disp('Generate SODF <> LS ...')
        f   = ['BaseA.SODF2LS.DV', num2str(ct), '.mat'];
        pf  = fullfile(InputData.PNAME_SOLUTION, f);
    end
    disp('file name is ... ')
    disp(f)
    
    A0  = sparse(size(q_DV,2),6*nnode);
    matlabpool(5)
    % STRAIN = {e11 e22 e33 sqrt(2)e12 sqrt(2)e13 sqrt(2)e23}
    % STRESS = {s11 s22 s33 sqrt(2)s12 sqrt(2)s13 sqrt(2)s23}
    if InputData.UseLSDF
        parfor k = 1:1:size(q_DV,2)
            if ~rem(k,500)
                disp([num2str(k/size(q_DV,2)*100) '%'])
            end
            
            A0(k,:) = BuildLsdfSpfMatrix(hkl_DV(:,k), ...
                frmesh, frmesh.symmetries, ...
                q_DV(:,k), InputData.fib_div, odf_data, ...
                w, InputData.deta, InputData.domega);
        end
    else
        parfor k = 1:1:size(q_DV,2)
            if ~rem(k,500)
                disp([num2str(k/size(q_DV,2)*100) '%'])
            end
            
            A0(k,:) = BuildSodfSpfMatrix(hkl_DV(:,k), ...
                frmesh, frmesh.symmetries, ...
                q_DV(:,k), InputData.fib_div, odf_data, ...
                w, InputData.S, InputData.deta, InputData.domega);
        end
    end
    matlabpool close
    disp(['saving ', f])
    save(pf, 'A0');
    toc
    
    ct  = ct + 1;
end