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

% XRDIMAGE.Image.samplename   = 'mishra_hr1_s1';
% FROOT   = {'mishra_hr1_s1_load11_pd'; ...
%     'mishra_hr1_s1_load12_pd'; ...
%     'mishra_hr1_s1_load13_pd'; ...
%     'mishra_hr1_s1_load14_pd'; ...
%     'mishra_hr1_s1_load15_pd'; ...
%     'mishra_hr1_s1_load16_pd' ; ...
%     'mishra_hr1_s1_load17_pd' ; ...
%     'mishra_hr1_s1_load18_pd' ; ...
%     'mishra_hr1_s1_load19_pd' ; ...
%     'mishra_hr1_s1_load20_pd' ; ...
%     'mishra_hr1_s1_load21_pd' ; ...
%     'mishra_hr1_s1_load22_pd' ; ...
%     'mishra_hr1_s1_load23_pd' ; ...
%     'mishra_hr1_s1_load24_pd' ; ...
%     'mishra_hr1_s1_load25_pd' ; ...
%     'mishra_hr1_s1_load26_pd' ; ...
%     'mishra_hr1_s1_load27_pd' ; ...
%     'mishra_hr1_s1_load28_pd' ; ...
%     'mishra_hr1_s1_load29_pd' ; ...
%     'mishra_hr1_s1_load30_pd' ; ...
%     'mishra_hr1_s1_load31_pd' ; ...
%     'mishra_hr1_s1_load32_pd' ; ...
%     };

% XRDIMAGE.Image.samplename   = 'mishra_hr2_s1';
% FROOT   = {'mishra_hr2_s1_load16_pd'; ...
%     'mishra_hr2_s1_load17_pd'; ...
%     'mishra_hr2_s1_load18_pd'; ...
%     'mishra_hr2_s1_load19_pd'; ...
%     'mishra_hr2_s1_load20_pd'; ...
%     'mishra_hr2_s1_load21_pd' ; ...
%     'mishra_hr2_s1_load22_pd' ; ...
%     'mishra_hr2_s1_load23_pd' ; ...
%     'mishra_hr2_s1_load24_pd' ; ...
%     'mishra_hr2_s1_load25_pd' ; ...
%     'mishra_hr2_s1_load26_pd' ; ...
%     'mishra_hr2_s1_load27_pd' ; ...
%     'mishra_hr2_s1_load28_pd' ; ...
%     'mishra_hr2_s1_load29_pd' ; ...
%     'mishra_hr2_s1_load30_pd' ; ...
%     'mishra_hr2_s1_load31_pd' ; ...
%     'mishra_hr2_s1_load32_pd' ; ...
%     'mishra_hr2_s1_load33_pd' ; ...
%     'mishra_hr2_s1_load34_pd' ; ...
%     'mishra_hr2_s1_load35_pd' ; ...
%     'mishra_hr2_s1_load36_pd' ; ...
%     'mishra_hr2_s1_load37_pd' ; ...
%     'mishra_hr2_s1_load38_pd' ; ...
%     'mishra_hr2_s1_load39_pd' ; ...
%     'mishra_hr2_s1_load40_pd' ; ...
%     'mishra_hr2_s1_load41_pd' ; ...
%     };

XRDIMAGE.Image.samplename   = 'mishra_hr3_s1';
FROOT   = {'mishra_hr3_s1_load16_pd'; ...
    'mishra_hr3_s1_load17_pd'; ...
    'mishra_hr3_s1_load18_pd'; ...
    'mishra_hr3_s1_load19_pd'; ...
    'mishra_hr3_s1_load20_pd'; ...
    'mishra_hr3_s1_load21_pd' ; ...
    'mishra_hr3_s1_load22_pd' ; ...
    'mishra_hr3_s1_load23_pd' ; ...
    'mishra_hr3_s1_load24_pd' ; ...
    'mishra_hr3_s1_load25_pd' ; ...
    'mishra_hr3_s1_load26_pd' ; ...
    'mishra_hr3_s1_load27_pd' ; ...
    'mishra_hr3_s1_load28_pd' ; ...
    'mishra_hr3_s1_load29_pd' ; ...
    'mishra_hr3_s1_load30_pd' ; ...
    'mishra_hr3_s1_load31_pd' ; ...
    'mishra_hr3_s1_load32_pd' ; ...
    'mishra_hr3_s1_load33_pd' ; ...
    'mishra_hr3_s1_load34_pd' ; ...
    'mishra_hr3_s1_load35_pd' ; ...
    'mishra_hr3_s1_load36_pd' ; ...
    };

%%% DATA REDUCTION FLAGS
Analysis_Options.make_polimg    = 1;
Analysis_Options.save_polimg    = 1;
Analysis_Options.fits_spectra   = 0;
Analysis_Options.save_fits      = 0;
Analysis_Options.find_instrpars = 0;
Analysis_Options.save_instrpars = 0;
Analysis_Options.find_detpars	= 0;
Analysis_Options.generateESG    = 1;

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
    
    %%% NEED TO LOOK AT THIS LAER
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
                case 'prrot'
                    warning('we have nothing here; maybe implement a flag that says prrot is angle or not?')
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
                                0, 0, XRDIMAGE.CakePrms)
                        case {'cor', 'cor32'}
                            %%% WritePolimg(pfname, polimg, instr, omega, chi, cakeprms)
                            WritePolimg(pfname_polimage, polimg, Instr, ...
                                omega_grid(jjj), chi_grid(jjj), XRDIMAGE.CakePrms)
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
                                0, 0, XRDIMAGE.CakePrms)
                        case {'cor', 'cor32'}
                            
                            %%% BuildESG(pfname, polImg, distance, wavelength, omega, chi, cakeParms)
                            BuildESG(pfname_esg, polimg, Instr.distance, Instr.wavelength, ...
                                omega_grid(jjj), chi_grid(jjj), XRDIMAGE.CakePrms)
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