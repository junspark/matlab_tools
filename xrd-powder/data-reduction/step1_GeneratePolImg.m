clear all
close all
clc

addpath(genpath('/home/beams/PARKJS/matlab/matlab_tools'));
% addpath(genpath('C:\Users\parkjs\Documents\GitHub\matlab_tools'));

%%% CALIBRATION FILES
XRDIMAGE.Calib.pname        = '/home/beams/S1IDUSER/mnt/orthros2/mpe1_oct20_bc';
XRDIMAGE.Calib.froot        = 'LaB6_1s';
XRDIMAGE.Calib.fnumber      = 6:7; % 4116 / 4117
XRDIMAGE.Calib.numdigs      = 6;
XRDIMAGE.Calib.fext         = 'ge3';
XRDIMAGE.Calib.corrected    = 'sum';

%%% SAMPLE FILES
XRDIMAGE.Image.pname        = '/home/beams/S1IDUSER/mnt/orthros2/mpe1_oct20_bc';
XRDIMAGE.Image.numdigs      = 6;
XRDIMAGE.Image.fext         = 'ge3';        %%% CHECK METADATA SEARCH
XRDIMAGE.Image.corrected    = 'cor32';
XRDIMAGE.Image.IsHydra      = 0;    % 0 = Single panel; 1 = GE1; 2 = GE2; 3 = GE3; 4 = GE4;

%%% METADATA FILES
XRDIMAGE.Image.per_scan_metadata    = '/home/beams/S1IDUSER/new_data/mpe1_oct20/saxs_waxs_fmt_fastpar.par';
XRDIMAGE.Image.per_frame_metadata   = '/home/beams/S1IDUSER/new_data/mpe1_oct20/mpe1_oct20_FF.par';

XRDIMAGE.Image.samplename   = 'hassani_sam1';
FROOT   = {'hassani_sam1_load0'; ...
    'hassani_sam1_load02_pd'; ...
    'hassani_sam1_load04_pd'; ...
    'hassani_sam1_load06_pd'; ...
    'hassani_sam1_load08_pd'; ...
    'hassani_sam1_load10_pd' ; ...
    'hassani_sam1_loadc1_pd' ; ...
    'hassani_sam1_loadc2_pd' ; ...
    'hassani_sam1_loadc3_pd' ; ...
    'hassani_sam1_loadc4_pd' ; ...
    'hassani_sam1_loadc5_pd' ; ...
    'hassani_sam1_loadc6_pd' ; ...
    'hassani_sam1_loadc7_pd' ; ...
    'hassani_sam1_loadc8_pd' ; ...
    'hassani_sam1_loadc9_pd' ; ...
    'hassani_sam1_loadc10_pd' ; ...
    'hassani_sam1_loadc11_pd' ; ...
    'hassani_sam1_loadc12_pd' ; ...
    'hassani_sam1_loadc13_pd' ; ...
    'hassani_sam1_loadc14_pd' ; ...
    'hassani_sam1_loadc15_pd' ; ...
    'hassani_sam1_loadc16_pd' ; ...
    'hassani_sam1_loadc17_pd' ; ...
    'hassani_sam1_loadc18_pd' ; ...
    'hassani_sam1_loadc19_pd' ; ...
    'hassani_sam1_loadc20_pd' ; ...
    'hassani_sam1_loadc21_pd' ; ...
    'hassani_sam1_loadc22_pd' ; ...
    'hassani_sam1_loadc23_pd' ; ...
    'hassani_sam1_loadc24_pd' ; ...
    'hassani_sam1_loadc25_pd' ; ...
    'hassani_sam1_loadc26_pd' ; ...
    'hassani_sam1_loadc27_pd' ; ...
    'hassani_sam1_loadc28_pd' ; ...
    'hassani_sam1_loadc29_pd' ; ...
    };

% XRDIMAGE.Image.samplename   = 'hassani_sam2';
% FROOT   = {'hassani_sam2_load00_pd'; ...
%     'hassani_sam2_load01_pd'; ...
%     'hassani_sam2_load02_pd'; ...
%     'hassani_sam2_load03_pd'; ...
%     'hassani_sam2_loadc01_pd'; ...
%     'hassani_sam2_loadc02_pd'; ...
%     'hassani_sam2_loadc03_pd'; ...
%     'hassani_sam2_loadc04_pd'; ...
%     'hassani_sam2_loadc05_pd'; ...
%     'hassani_sam2_loadc06_pd'; ...
%     'hassani_sam2_loadc07_pd'; ...
%     'hassani_sam2_loadc08_pd'; ...
%     'hassani_sam2_loadc09_pd'; ...
%     };

for iiii = 1:1:length(FROOT)
    %%% INPUT PARAMETERS
    XRDIMAGE.Image.froot        = FROOT{iiii};     %%% FIND FNUM BASED ON THIS FILE NAME ROOT 32
    
    %%% GET IMAGE NUMBERS
    metadata_per_scan   = ReadSpecParFile(XRDIMAGE.Image.per_scan_metadata, 'Version', 'saxs_waxs_fmt_fastpar_v3');
    metadata_per_frame  = ReadSpecParFile(XRDIMAGE.Image.per_frame_metadata, 'Version' ,'mpe_ff_per_frame_v2');
    idx_froot           = ismember(metadata_per_scan.det3_fname, XRDIMAGE.Image.froot);
    
    XRDIMAGE.Image.fnumber      = metadata_per_scan.det3_fnum(idx_froot) - 1;
    switch XRDIMAGE.Image.corrected
        case 'sum'
            XRDIMAGE.Image.scan_nframes = 1;
        case {'cor', 'cor32'}
            XRDIMAGE.Image.scan_nframes = metadata_per_scan.scan_nframes(idx_froot);
            XRDIMAGE.Image.scan_mtr     = metadata_per_scan.scan_mtr(idx_froot);
            XRDIMAGE.Image.scan_ini     = metadata_per_scan.scan_ini(idx_froot);
            XRDIMAGE.Image.scan_fin     = metadata_per_scan.scan_fin(idx_froot);
        otherwise
            warning('num frames per scan could not be determined')
    end
    
    numimg  = length(XRDIMAGE.Calib.fnumber);
    for iii = 1:1:numimg
        fname   = sprintf('%s_%06d.%s.%s', ...
            XRDIMAGE.Calib.froot, XRDIMAGE.Calib.fnumber(iii), XRDIMAGE.Calib.fext, XRDIMAGE.Calib.corrected);
        pfname{iii, 1}  = fullfile(XRDIMAGE.Calib.pname, ...
            XRDIMAGE.Calib.froot, XRDIMAGE.Calib.fext, fname);
    end
    
    %%% INSTRUMENT PARAMETERS GETS LOADED
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
    
    Instr.omega = 0;
    Instr.chi   = 0;
    
    XRDIMAGE.Instr  = Instr;
    
    %%% CAKE PARAMETERS
    XRDIMAGE.CakePrms.bins(1)   = 36;           % number of azimuthal bins
    XRDIMAGE.CakePrms.bins(2)   = 1000;         % number of radial bins
    XRDIMAGE.CakePrms.bins(3)   = 70;           % number of angular bins
    
    XRDIMAGE.CakePrms.origin(1) = 1024.190;     % apparent X center in pixels // THIS IS WHAT YOU SEE ON FIGURE 1
    XRDIMAGE.CakePrms.origin(2) = 1062.110;     % apparent Y center in pixels // THIS IS WHAT YOU SEE ON FIGURE 1
    XRDIMAGE.CakePrms.origin(2) = XRDIMAGE.Instr.numpixelsVert-XRDIMAGE.CakePrms.origin(2); %%% CONVERT TO IMAGE COORDIANTES
    
    XRDIMAGE.CakePrms.sector(1) = -360/XRDIMAGE.CakePrms.bins(1)/2;     % start azimuth (min edge of bin) in degrees
    XRDIMAGE.CakePrms.sector(2) = 360-360/XRDIMAGE.CakePrms.bins(1)/2;  % stop  azimuth (max edge of bin) in degrees
    XRDIMAGE.CakePrms.sector(3) = 150;  % start radius (min edge of bin) in pixels
    XRDIMAGE.CakePrms.sector(4) = 1000;  % stop  radius (max edge of bin) in pixels
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eta_step    = (XRDIMAGE.CakePrms.sector(2) - XRDIMAGE.CakePrms.sector(1))/XRDIMAGE.CakePrms.bins(1);
    eta_ini     = XRDIMAGE.CakePrms.sector(1) + eta_step/2;
    eta_fin     = XRDIMAGE.CakePrms.sector(2) - eta_step/2;
    azim        = eta_ini:eta_step:eta_fin;
    XRDIMAGE.CakePrms.azim      = azim;
    XRDIMAGE.CakePrms.fastint   = 1;
    
    %%% DATA REDUCTION FLAGS
    Analysis_Options.make_polimg    = 1;
    Analysis_Options.save_polimg    = 1;
    Analysis_Options.fits_spectra   = 0;
    Analysis_Options.save_fits      = 0;
    Analysis_Options.find_instrpars = 0;
    Analysis_Options.save_instrpars = 0;
    Analysis_Options.find_detpars	= 0;
    Analysis_Options.generateESG    = 1;
    
    %%% PK FITTING OPTIONS
    Analysis_Options.PkFitOptions   = optimset(...
        'MaxIter', 5e5,...
        'MaxFunEvals',3e5);
    
    Analysis_Options.InstrPrmFitOptions = optimset(...
        'DerivativeCheck', 'off', ...
        'MaxIter', 1e5, ...
        'MaxFunEvals', 3e5, ...
        'TypicalX',[100 -100 1000 0.1 0.1 XRDIMAGE.Instr.detpars], ...
        'Display','final');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% GENERATE MESH FOR INTEGRATION
    %%% IF POLIMG NEEDS TO BE GENERATED
    if Analysis_Options.make_polimg
        if ~XRDIMAGE.CakePrms.fastint
            DetectorMesh    = BuildMeshDetector(XRDIMAGE.Instr.numpixels, XRDIMAGE.CakePrms);
        else
            DetectorMesh    = 0;
        end
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
                for jjj = 1:1:XRDIMAGE.Image.scan_nframes
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
    %%% LOAD XRD IMAGES & GENERATE POLIMG FILE
    if Analysis_Options.make_polimg
        delete(gcp)
        parpool(10)
        parfor iii = 1:1:numimg
            scan_mtr    = XRDIMAGE.Image.scan_mtr{iii};
            switch scan_mtr
                case 'aero'
                    omega_grid  = linspace(XRDIMAGE.Image.scan_ini(iii), ...
                        XRDIMAGE.Image.scan_fin(iii), ...
                        XRDIMAGE.Image.scan_nframes(iii)+1);
                    omega_grid  = (omega_grid(1:end-1) + omega_grid(2:end))./2;
                    
                    chi_grid    = 0.*omega_grid;
                    prrot_grid  = 0.*omega_grid;
                case 'prrot'
                    prrot_grid  = linspace(XRDIMAGE.Image.scan_ini(iii), ...
                        XRDIMAGE.Image.scan_fin(iii), ...
                        XRDIMAGE.Image.scan_nframes(iii)+1);
                    prrot_grid  = (prrot_grid(1:end-1) + prrot_grid(2:end))./2;
                    
                    omega_grid  = 0.*prrot_grid;
                    chi_grid    = 0.*prrot_grid;
            end
            
            for jjj = 1:1:XRDIMAGE.Image.scan_nframes
                tic
                disp('###########################')
                disp(sprintf('Looking at %s', pfname{iii, jjj}))
                disp('###########################')
                
                pfname_polimage = [pfname{iii,jjj}, '.polimg.mat'];
                pfname_esg      = [pfname{iii,jjj}, '.esg'];
                
                imgi    = ReadSUM(pfname{iii,1});
                
                %         figure(1)
                %         hold off
                %         imagesc(rot90(imgi,1))
                %         caxis([-10 3000])
                %         axis equal tight
                %         colorbar vert
                %         hold on
                %         xlabel('X_L (pixels)')
                %         ylabel('Y_L (pixels)')
                %         title('Ensure that image matches the coordinate system')
                %         text(0, 0, 'TO')
                %         text(2048, 0, 'TI')
                %         text(0, 0, 'TO')
                %         text(0, 2048, 'BO')
                %         text(2048, 2048, 'BI')
                
                %%% POLAR REBINNING
                polimg  = PolarBinXRD(DetectorMesh, ...
                    XRDIMAGE.Instr, ...
                    XRDIMAGE.CakePrms, ...
                    imgi, 'PlotProgress', 'off');
                
                mapped_tth  = GeometricModelXRDSwitch(XRDIMAGE.Instr, polimg);
                polimg.mapped_tth_for_intensity = mapped_tth;
                
                [tth_grid, d_grid, intensity_in_tth_grid]   = MapIntensityToTThGrid(XRDIMAGE, polimg);
                polimg.tth_grid                 = tth_grid;
                polimg.d_grid                   = d_grid;
                polimg.intensity_in_tth_grid    = intensity_in_tth_grid;
                
                if Analysis_Options.save_polimg
                    disp(sprintf('Saving polimg for %s', pfname{iii,jjj}))
                    % save(pfname_polimage, 'polimg', 'XRDIMAGE')
                    
                    switch XRDIMAGE.Image.corrected
                        case 'sum'
                            %%% WritePolimg(pfname, polimg, instr, omega, chi, cakeprms)
                            WritePolimg(pfname_polimage, polimg, Instr, ...
                                0, 0, 0, XRDIMAGE.CakePrms)
                        case {'cor', 'cor32'}
                            %%% WritePolimg(pfname, polimg, instr, omega, chi, cakeprms)
                            WritePolimg(pfname_polimage, polimg, Instr, ...
                                omega_grid(jjj), chi_grid(jjj), prrot_grid(jjj), XRDIMAGE.CakePrms)
                    end
                else
                    disp(sprintf('Not saving polimg for %s', pfname{iii,jjj}))
                end
                
                if Analysis_Options.generateESG
                    disp(sprintf('Saving MAUD compatible ESG file for %s', pfname{iii,jjj}))
                    switch XRDIMAGE.Image.corrected
                        case 'sum'
                            %%% BuildESG(pfname, polImg, distance, wavelength, omega, chi, cakeParms)
                            BuildESG(pfname_esg, polimg, Instr.distance, Instr.wavelength, ...
                                0, 0, 0, XRDIMAGE.CakePrms)
                        case {'cor', 'cor32'}
                            %%% BuildESG(pfname, polImg, distance, wavelength, omega, chi, cakeParms)
                            BuildESG(pfname_esg, polimg, Instr.distance, Instr.wavelength, ...
                                omega_grid(jjj), chi_grid(jjj), prrot_grid(jjj), XRDIMAGE.CakePrms)
                    end
                else
                    disp(sprintf('Not saving esg for %s', pfname{iii,jjj}))
                end
                
                %         figure(2)
                %         subplot(1,2,1)
                %         imagesc(log(polimg.intensity)), axis square tight
                %         hold off
                %
                %         subplot(1,2,2)
                %         plot(polimg.radius, polimg.intensity)
                %         hold off
                %         disp(' ')
                %
                %         figure(3)
                %         imagesc(log(abs(polimg.intensity_in_tth_grid))), axis square tight
                %         title('Caked image // radial position is corrected')
                %         hold off
                toc
            end
        end
    end
end