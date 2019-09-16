clear all
close all
clc

%%% INPUT PARAMETERS
XRDIMAGE.Image.pname        = '/home/beams/S1IDUSER/mnt/s1b/__eval/projects_parkjs/dexela_workflow/shastri_jul19/dunand_AlCe_alloy/data/';
XRDIMAGE.Image.fbase        = 'AlCe_sam2_load40';
XRDIMAGE.Image.fnumber      = 1876:2175; % 4116 / 4117
XRDIMAGE.Image.numframe     = 1;
XRDIMAGE.Image.numdigs      = 6;
XRDIMAGE.Image.fext         = 'tif';
XRDIMAGE.Image.corrected    = 0;
XRDIMAGE.Image.IsHydra      = 0;    % 0 = Single panel; 1 = GE1; 2 = GE2; 3 = GE3; 4 = GE4;

%%% DARK FILES ONLY USED IF THE IMAGES ARE UNCORRECTED
XRDIMAGE.DarkField.pname    = fullfile(XRDIMAGE.Image.pname, XRDIMAGE.Image.fbase);
XRDIMAGE.DarkField.fbase    = 'dark_before';
XRDIMAGE.DarkField.fnumber  = 1866:1875;
XRDIMAGE.DarkField.numframe = 1;
XRDIMAGE.DarkField.numdigs  = 6;
XRDIMAGE.DarkField.fext     = 'tif';

%%% CALIBRATION FILE
XRDIMAGE.Calib.pname        = '/home/beams/S1IDUSER/mnt/s1b/__eval/projects_parkjs/dexela_workflow/shastri_jul19/dunand_AlCe_alloy/data/';
XRDIMAGE.Calib.fbase        = 'CeO2_0pt5s';
XRDIMAGE.Calib.fnumber      = 53;
XRDIMAGE.Calib.fext         = 'tif';

%%% METADATA FILE
XRDIMAGE.MetaData.pfname    = '/home/beams/S1IDUSER/new_data/shastri_jul19/saxs_waxs_fmt_fastpar2.par';
XRDIMAGE.MetaData.data      = ReadSpecParFile(XRDIMAGE.MetaData.pfname, ...
    'Version', 'saxs_waxs_fmt_fastpar_v1');

fname_pattern   = sprintf('%%s_%%0%dd.%%s', XRDIMAGE.Image.numdigs);

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
for i = 1:1:length(XRDIMAGE.Calib.fnumber)
    fname_instr     = [sprintf(fname_pattern, XRDIMAGE.Calib.fbase, XRDIMAGE.Calib.fnumber(i), XRDIMAGE.Calib.fext), '.instr.mat'];
    pfname_instr	= fullfile(XRDIMAGE.Calib.pname, XRDIMAGE.Calib.fbase, fname_instr);
    
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
XRDIMAGE.CakePrms.bins(1)   = 72;           % number of azimuthal bins
XRDIMAGE.CakePrms.bins(2)   = 3500;         % number of radial bins
XRDIMAGE.CakePrms.bins(3)   = 60;            % number of angular bins

XRDIMAGE.CakePrms.origin(1) = 1942.795;   % apparent X center in pixels // THIS IS WHAT YOU SEE ON FIGURE 1
XRDIMAGE.CakePrms.origin(2) = 1555.548;     % apparent Y center in pixels // THIS IS WHAT YOU SEE ON FIGURE 1

XRDIMAGE.CakePrms.sector(1) = -360/XRDIMAGE.CakePrms.bins(1)/2;     % start azimuth (min edge of bin) in degrees
XRDIMAGE.CakePrms.sector(2) = 360-360/XRDIMAGE.CakePrms.bins(1)/2;  % stop  azimuth (max edge of bin) in degrees
XRDIMAGE.CakePrms.sector(3) = 150;      % start radius (min edge of bin) in pixels
XRDIMAGE.CakePrms.sector(4) = 1150;     % stop  radius (max edge of bin) in pixels

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
%%% LOAD XRD IMAGES
%%% BACKGROUND
if XRDIMAGE.Image.corrected
    disp('###########################')
    fprintf('images are already corrected for background.\n');
    disp('###########################')
else
    disp('###########################')
    fprintf('loading background file for dark.\n');
    disp('###########################')
    
    for i = 1:1:length(XRDIMAGE.DarkField.fnumber)
        fname	= sprintf(fname_pattern, XRDIMAGE.DarkField.fbase, XRDIMAGE.DarkField.fnumber(i), XRDIMAGE.DarkField.fext);
        pfname  = fullfile(XRDIMAGE.DarkField.pname, fname);
        
        if i == 1
            bkg = ReadDexela(pfname, 'SPECSettingCorrect', false, 'Orientation', 'cable_UP');
        else
            bkg = bkg + ReadDexela(pfname, 'SPECSettingCorrect', false, 'Orientation', 'cable_UP');
        end
    end
    bkg = bkg./length(XRDIMAGE.DarkField.fnumber);
        
    numpixHorz  = size(bkg,1);
    numpixVert  = size(bkg,2);
    
    XRDIMAGE.Instr.numpixelsHorz    = numpixHorz;
    XRDIMAGE.Instr.numpixelsVert    = numpixVert;
    
    if numpixHorz == 3888
        XRDIMAGE.Instr.pixelsizeHorz    = 0.0748;          % mm
    elseif numpixHorz == 1944
        XRDIMAGE.Instr.pixelsizeHorz    = 0.1496;          % mm
    elseif numpixHorz == 972
        XRDIMAGE.Instr.pixelsizeHorz    = 0.2992;          % mm
    end
    
    if numpixVert == 3072
        XRDIMAGE.Instr.pixelsizeVert    = 0.0748;          % mm
    elseif numpixVert == 1536
        XRDIMAGE.Instr.pixelsizeVert    = 0.1496;          % mm
    elseif numpixVert == 768
        XRDIMAGE.Instr.pixelsizeVert    = 0.2992;          % mm
    end
    
    XRDIMAGE.Instr.detsizeHorz  = XRDIMAGE.Instr.pixelsizeHorz*XRDIMAGE.Instr.numpixelsHorz;        % mm
    XRDIMAGE.Instr.detsizeVert  = XRDIMAGE.Instr.pixelsizeVert*XRDIMAGE.Instr.numpixelsVert;        % mm
    
    %%% CONVERT TO IMAGE COORDINATES
    XRDIMAGE.CakePrms.origin(2) = XRDIMAGE.Instr.numpixelsVert - XRDIMAGE.CakePrms.origin(2);
end

if Analysis_Options.make_polimg
    if ~XRDIMAGE.CakePrms.fastint
        DetectorMesh    = BuildMeshDetector(XRDIMAGE.Instr.numpixelsHorz, XRDIMAGE.Instr.numpixelsVert, XRDIMAGE.CakePrms);
    else
        DetectorMesh    = 0;
    end
end

numimg  = length(XRDIMAGE.Image.fnumber);
if Analysis_Options.make_polimg
    for i = 1:1:numimg
        fname   = sprintf(fname_pattern, XRDIMAGE.Image.fbase, XRDIMAGE.Image.fnumber(i), XRDIMAGE.Image.fext);
        pfname  = fullfile(XRDIMAGE.Image.pname, XRDIMAGE.Image.fbase, fname);
        disp('###########################')
        disp(sprintf('Looking at %s', pfname))
        disp('###########################')
        
        pfname_polimage = [pfname, '.polimg.mat'];
        pfname_esg      = [pfname, '.esg'];
        
        if XRDIMAGE.Image.corrected
            imgi    = ReadSUM(pfname{i,1});
        else
            imgi    = bkg.*0;
            for j = 1:1:XRDIMAGE.Image.numframe
                imgj    = ReadDexela(pfname, 'SPECSettingCorrect', false, 'Orientation', 'cable_UP');
                imgi    = imgi + imgj;
            end
            imgi    = imgi - bkg.*XRDIMAGE.Image.numframe;
        end
        
        %%% FIND THE IMAGE STACK IN THE METADATA
        idx = find(XRDIMAGE.Image.fnumber(i) >= XRDIMAGE.MetaData.data.imnum_ini);
        idx = idx(end);
        if ((XRDIMAGE.Image.fnumber(i) >= XRDIMAGE.MetaData.data.imnum_ini(idx)) && ...
            (XRDIMAGE.Image.fnumber(i) <= XRDIMAGE.MetaData.data.imnum_fin(idx)))
            disp('found the image in the metadata file. continue with polimg generation.')
            disp('update XRDIMAGE.Instr with angle information')
            
            omeini  = double(XRDIMAGE.MetaData.data.scan_ini(idx));
            omefin  = double(XRDIMAGE.MetaData.data.scan_fin(idx));
            nframes = double(XRDIMAGE.MetaData.data.scan_nframes(idx));
            omestep = (omefin - omeini)/nframes;
            
            omega   = omeini + omestep*(XRDIMAGE.Image.fnumber(i) - double(XRDIMAGE.MetaData.data.imnum_ini(idx)));
            chi     = 0;
            
            XRDIMAGE.Instr.omega    = omega;
            XRDIMAGE.Instr.chi      = chi;
        else
            disp('did not find the image in the metadata file.')
            return
        end
        
        %%% SYNTHETIC IMAGE
%         xxx = (1:1:3888)./388.8;
%         yyy = (1:1:3072)./(307.2*2);
%         imgi    = xxx'*yyy;
%         imgi(2000:2200, :) = 27;
        
        figure(1)
        hold off
        imagesc(rot90(imgi, 1)) %% DO NOT CHANGE
        caxis([-10 3000])
        axis equal tight
        colorbar vert
        hold on
        xlabel('X_L (pixels)')
        ylabel('Y_L (pixels)')
        title('Ensure that image matches the coordinate system')
        text(XRDIMAGE.Instr.numpixelsHorz, 0, 'TI')
        text(0, 0, 'TO')
        text(0, XRDIMAGE.Instr.numpixelsVert, 'BO')
        text(XRDIMAGE.Instr.numpixelsHorz, XRDIMAGE.Instr.numpixelsVert, 'BI')
        
        %%% POLAR REBINNING
        polimg	= PolarBinXRD(DetectorMesh, ...
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
            disp(sprintf('Saving polimg for %s', pfname))
            save(pfname_polimage, 'polimg', 'XRDIMAGE')
        else
            disp(sprintf('Not saving polimg for %s', pfname))
        end
        
        if Analysis_Options.generateESG
            disp(sprintf('Saving MAUD compatible ESG file for %s', pfname))
            BuildESG(pfname_esg, polimg, ...
                XRDIMAGE.Instr.distance, ...
                XRDIMAGE.Instr.wavelength, ...
                XRDIMAGE.Instr.omega, ...
                XRDIMAGE.Instr.chi, XRDIMAGE.CakePrms)
        end
        
        figure(2)
        subplot(1,2,1)
        imagesc(log(abs(polimg.intensity))), axis square tight
        hold off
        
        subplot(1,2,2)
        plot(polimg.radius, polimg.intensity)
        hold off
        disp(' ')
    end
end