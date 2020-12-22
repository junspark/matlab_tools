clear all
close all
clc

XRDIMAGE.ExpID              = 'lywang_nov20';

% addpath(genpath('/home/beams/PARKJS/matlab/matlab_tools'));
addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

%%% CALIBRATION FILES
XRDIMAGE.Calib.pname        = 'C:\Users\parkjs\Documents\MATLAB\lywang_nov20\polimg';
XRDIMAGE.Calib.froot        = 'LaB6_1s_Mg';
XRDIMAGE.Calib.fnumber      = 8:9; % 4116 / 4117
XRDIMAGE.Calib.numdigs      = 6;
XRDIMAGE.Calib.fext         = 'ge3';
XRDIMAGE.Calib.corrected    = 'sum';

%%% SAMPLE FILES
XRDIMAGE.Image.pname        = 'C:\Users\parkjs\Documents\MATLAB\lywang_nov20\polimg';
XRDIMAGE.Image.numdigs      = 6;
XRDIMAGE.Image.fext         = 'ge3';        %%% CHECK METADATA SEARCH
XRDIMAGE.Image.corrected    = 'cor32';
XRDIMAGE.Image.IsHydra      = 0;    % 0 = Single panel; 1 = GE1; 2 = GE2; 3 = GE3; 4 = GE4;

%%% METADATA FILES
XRDIMAGE.Image.per_scan_metadata        = 'C:\Users\parkjs\Documents\MATLAB\lywang_nov20\lywang_nov20_metadata\saxs_waxs_fmt_fastpar.par';
XRDIMAGE.Image.per_scan_waxs_metadata   = 'C:\Users\parkjs\Documents\MATLAB\lywang_nov20\lywang_nov20_metadata\waxs_hydra_exposures.par';
XRDIMAGE.Image.per_frame_metadata       = 'C:\Users\parkjs\Documents\MATLAB\lywang_nov20\lywang_nov20_metadata\lywang_nov20_FF.par';
XRDIMAGE.Image.per_frame_waxs_metadata  = 'C:\Users\parkjs\Documents\MATLAB\lywang_nov20\lywang_nov20_metadata\per_frame_waxs_hydra_exposures2.par';

numimg  = length(XRDIMAGE.Calib.fnumber);
for iii = 1:1:numimg
    fname   = sprintf('%s_%06d.%s.%s', ...
        XRDIMAGE.Calib.froot, XRDIMAGE.Calib.fnumber(iii), XRDIMAGE.Calib.fext, XRDIMAGE.Calib.corrected);
    pfname{iii, 1}  = fullfile(XRDIMAGE.Calib.pname, ...
        XRDIMAGE.Calib.froot, XRDIMAGE.Calib.fext, fname);
end

%%% INSTRUMENT PARAMETERS GETS LOADED TO DO SOME INITIAL CALCULATION
Instr.energy        = 0;
Instr.wavelength	= 0;
Instr.distance      = 0;
Instr.centers       = [0, 0];
Instr.gammaX        = 0;
Instr.gammaY        = 0;
Instr.detsizeHorz	= 0;
Instr.detsizeVert   = 0;
Instr.pixelsizeHorz	= 0;
Instr.pixelsizeVert	= 0;
Instr.numpixelsHorz	= 0;
Instr.numpixelsVert	= 0;
Instr.imrotation    = 0;
Instr.dettype       = '2a';
Instr.detpars       = [0 0 0 0 0 0];
for iii = 1:1:length(XRDIMAGE.Calib.fnumber)
    pfname_instr    = [pfname{iii,1}, '.instr.mat'];
    Instri  = load(pfname_instr);
    Instri  = Instri.Instr;
    
    Instr.energy        = Instr.energy + Instri.energy;
    Instr.wavelength    = Instr.wavelength + Instri.wavelength;
    Instr.distance      = Instr.distance + Instri.distance;
    Instr.centers       = Instr.centers + Instri.centers;
    Instr.gammaX        = Instr.gammaX + Instri.gammaX;
    Instr.gammaY        = Instr.gammaY + Instri.gammaY;
    Instr.detsizeHorz   = Instr.detsizeHorz + Instri.detsizeHorz;
    Instr.detsizeVert   = Instr.detsizeVert + Instri.detsizeVert;
    Instr.pixelsizeHorz	= Instr.pixelsizeHorz + Instri.pixelsizeHorz;
    Instr.pixelsizeVert	= Instr.pixelsizeVert + Instri.pixelsizeVert;
    Instr.numpixelsHorz	= Instr.numpixelsHorz + Instri.numpixelsHorz;
    Instr.numpixelsVert	= Instr.numpixelsVert + Instri.numpixelsVert;
    Instr.detpars       = Instr.detpars + Instri.detpars;
    Instr.imrotation    = 0;
    Instr.dettype       = '2a';
end
Instr.energy        = Instr.energy./length(XRDIMAGE.Calib.fnumber);
Instr.wavelength    = Instr.wavelength./length(XRDIMAGE.Calib.fnumber);
Instr.distance      = Instr.distance./length(XRDIMAGE.Calib.fnumber);
Instr.centers       = Instr.centers./length(XRDIMAGE.Calib.fnumber);
Instr.gammaX        = Instr.gammaX./length(XRDIMAGE.Calib.fnumber);
Instr.gammaY        = Instr.gammaY./length(XRDIMAGE.Calib.fnumber);
Instr.detsizeHorz   = Instr.detsizeHorz./length(XRDIMAGE.Calib.fnumber);
Instr.detsizeVert   = Instr.detsizeVert./length(XRDIMAGE.Calib.fnumber);
Instr.pixelsizeHorz	= Instr.pixelsizeHorz./length(XRDIMAGE.Calib.fnumber);
Instr.pixelsizeVert	= Instr.pixelsizeVert./length(XRDIMAGE.Calib.fnumber);
Instr.numpixelsHorz	= Instr.numpixelsHorz./length(XRDIMAGE.Calib.fnumber);
Instr.numpixelsVert	= Instr.numpixelsVert./length(XRDIMAGE.Calib.fnumber);
Instr.detpars       = Instr.detpars./length(XRDIMAGE.Calib.fnumber);

XRDIMAGE.Image.samplename   = 'Ti_EOS_H_Ti44_400_sam1_waxs_cont_loading';
FROOT   = {'Ti_EOS_H_Ti44_400_sam1_waxs_cont_loading'; ...
    };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATERIAL PARAMETERS - 
Material.num        = 1;
Material.lattparms  = [2.92 4.66];
Material.structure  = 'hcp';
Material.hkls       = load([Material.structure, '.hkls']);

%%% CALCULATE THEORETICAL TTH
[d, th] = PlaneSpacings(Material.lattparms, ...
    'hexagonal', Material.hkls', ...
    Instr.wavelength);
tth     = 2*th;
d_spacing_range = 0.03;
d_spacing_UB    = (1 + d_spacing_range)*d;
d_spacing_LB    = (1 - d_spacing_range)*d;

tth_UB  = 2.*asind(Instr.wavelength/2)./d_spacing_LB;
tth_LB  = 2.*asind(Instr.wavelength/2)./d_spacing_UB;

Material.tth        = tth;
Material.tth_UB     = tth_UB;
Material.tth_LB     = tth_LB;

Material.d_spacing  = d;
Material.pkidx      = {...
    [1] [2 3] [4]
    };

Material.numbounds  = length(Material.pkidx);
numpk   = 0;
for iii = 1:1:Material.numbounds
    Material.pkrange(:,iii)  = [ ...
        min(tth_LB(Material.pkidx{iii})); ...
        max(tth_UB(Material.pkidx{iii})); ...
        ];
    numpk   = numpk + length(Material.pkidx{iii});
end
Material.numpk  = numpk;
Material.pkbck  = 2;
Material.pkfunc = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PK FITTING OPTIONS
Analysis_Options.PkFuncOptions.pfunc_type	= 'pseudoVoigt';
Analysis_Options.PkFuncOptions.pbkg_order	= 2;
Analysis_Options.PkFitOptimizationOptions   = optimset(...
    'MaxIter', 5e5, ...
    'MaxFunEvals',3e5);

%%% DATA REDUCTION FLAGS
Analysis_Options.make_polimg            = 0;
Analysis_Options.save_polimg            = 0;
Analysis_Options.fits_spectra           = 1;
Analysis_Options.use_parallel_process   = 1;
Analysis_Options.save_fits              = 1;
Analysis_Options.find_instrpars         = 0;
Analysis_Options.save_instrpars         = 0;
Analysis_Options.find_detpars           = 0;
Analysis_Options.generateESG            = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% GET IMAGE NUMBERS
metadata_per_scan       = ReadSpecParFile(XRDIMAGE.Image.per_scan_metadata, 'Version', 'saxs_waxs_fmt_fastpar_v3');
metadata_per_scan_waxs  = ReadSpecParFile(XRDIMAGE.Image.per_scan_waxs_metadata, 'Version' ,'mpe_standard');
metadata_per_frame      = ReadSpecParFile(XRDIMAGE.Image.per_frame_metadata, 'Version' ,'mpe_ff_per_frame_v2');
metadata_per_frame_waxs = ReadSpecPerFrameParFile(XRDIMAGE.Image.per_frame_waxs_metadata);

for iiii = 1:1:length(FROOT)
    %%% INPUT PARAMETERS
    XRDIMAGE.Image.froot    = FROOT{iiii};     %%% FIND FNUM BASED ON THIS FILE NAME ROOT 32
    
    idx_froot_in_metadata_per_scan          = ismember(metadata_per_scan.det3_fname, XRDIMAGE.Image.froot);
    idx_froot_in_metadata_per_scan_waxs     = ismember(metadata_per_scan_waxs.det3_fname, XRDIMAGE.Image.froot);
    idx_froot_in_metadata_per_frame         = ismember(metadata_per_frame.det_fname, XRDIMAGE.Image.froot);
    idx_froot_in_metadata_per_frame_waxs    = ismember(metadata_per_frame_waxs.det_fname, XRDIMAGE.Image.froot);
    
    XRDIMAGE.Image.fnumber      = metadata_per_scan.det3_fnum(idx_froot_in_metadata_per_scan) - 1;
    
    switch XRDIMAGE.Image.corrected
        case 'sum'
            XRDIMAGE.Image.scan_nframes = 1;
        case {'cor', 'cor32'}
            XRDIMAGE.Image.scan_nframes = metadata_per_frame_waxs.det_frames_per_file(idx_froot_in_metadata_per_frame_waxs);
            
            XRDIMAGE.Image.scan_mtr     = metadata_per_scan.scan_mtr(idx_froot_in_metadata_per_scan);
            XRDIMAGE.Image.scan_ini     = metadata_per_scan.scan_ini(idx_froot_in_metadata_per_scan);
            XRDIMAGE.Image.scan_fin     = metadata_per_scan.scan_fin(idx_froot_in_metadata_per_scan);
        otherwise
            warning('num frames per scan could not be determined')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% GENERATE LIST OF FILES TO LOOK AT
    numimg  = length(XRDIMAGE.Image.fnumber);
    for iii = 1:1:numimg
        switch XRDIMAGE.Image.corrected
            case 'sum'
                fname   = sprintf('%s_%06d.%s.%s', ...
                    XRDIMAGE.Image.froot, XRDIMAGE.Image.fnumber(iii), XRDIMAGE.Image.fext, XRDIMAGE.Image.corrected);
                pfname{iii, 1}  = fullfile(XRDIMAGE.Image.pname, ...
                    XRDIMAGE.Image.samplename, XRDIMAGE.Image.fext, fname);
            case {'cor', 'cor32'}
%                 fname   = sprintf('%s_%06d.%s', ...
%                     XRDIMAGE.Image.froot, XRDIMAGE.Image.fnumber(iii), XRDIMAGE.Image.fext);
%                 pfname_ge   = fullfile('/home/beams/S1IDUSER/mnt/s1c/lywang_nov20/ge3/', fname);
%                 
%                 XRDIMAGE.Image.scan_nframes(iii)    = CalcNumFramesGE(pfname_ge);
                
                for jjj = 1:1:XRDIMAGE.Image.scan_nframes(iii)
                    fname   = sprintf('%s_%06d.%s_frame_%d.%s', ...
                        XRDIMAGE.Image.froot, XRDIMAGE.Image.fnumber(iii), XRDIMAGE.Image.fext, jjj, XRDIMAGE.Image.corrected);
                    pfname{iii, jjj}    = fullfile(XRDIMAGE.Image.pname, ...
                        XRDIMAGE.Image.samplename, XRDIMAGE.Image.fext, fname);
                end
            otherwise
                warning('check XRDIMAGE.Image.corrected paramter')
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% LOAD XRD POLIMG & FIT PEAKS
    if Analysis_Options.fits_spectra
        for iii = 1:1:numimg
            
            idx_fnum_per_frame_waxs = find(metadata_per_frame_waxs.det_fnum == XRDIMAGE.Image.fnumber(iii));
            
            % metadata_per_frame_waxs.per_frame_data0{idx_fnum_per_frame_waxs} %%% DISP (um)
            % metadata_per_frame_waxs.per_frame_data1{idx_fnum_per_frame_waxs} %%% I0
            % metadata_per_frame_waxs.per_frame_data2{idx_fnum_per_frame_waxs} %%% I1
            % metadata_per_frame_waxs.per_frame_data3{idx_fnum_per_frame_waxs}  %%% LOAD (N)
            per_frame_data0 = metadata_per_frame_waxs.per_frame_data0{idx_fnum_per_frame_waxs};
            per_frame_data3 = metadata_per_frame_waxs.per_frame_data3{idx_fnum_per_frame_waxs};
            
            MetaDataFieldName	= {'Load', 'Disp'};
            
            if license('test','Parallel_Computing_Toolbox') && Analysis_Options.use_parallel_process
                delete(gcp)
                if XRDIMAGE.Image.scan_nframes(iii) > 20
                    parpool(20);
                else
                    parpool(XRDIMAGE.Image.scan_nframes(iii));
                end
                
                parfor jjj = 1:1:XRDIMAGE.Image.scan_nframes(iii)
                    tic
                    FitPeaksPerPolImage(pfname{iii, jjj}, Material, Analysis_Options, 'MetaDataFieldName', MetaDataFieldName, ...
                        'MetaDataFieldValues', [per_frame_data0(jjj), per_frame_data3(jjj)]);
                    toc
                end
                
                delete(gcp);
            else
                fprintf('parallel computing toolbox does not exist ...\n')
                fprintf('continue with serial operation.\n')
                
                for jjj = 1:1:XRDIMAGE.Image.scan_nframes(iii)
                    tic
                    FitPeaksPerPolImage(pfname{iii, jjj}, Material, Analysis_Options, 'MetaDataFieldName', MetaDataFieldName, ...
                        'MetaDataFieldValues', [per_frame_data0(jjj), per_frame_data3(jjj)]);
                    toc
                end
            end
        end
    end
end