clear all
close all
clc

%%% INPUT PARAMETERS
XRDIMAGE.Image.pname        = '/home/beams/S1IDUSER/mnt/s1b/__eval/projects_parkjs/dexela_workflow/shastri_jul19/dunand_AlCe_alloy/data/CeO2_0pt5s';
XRDIMAGE.Image.fbase        = 'CeO2_0pt5s';
XRDIMAGE.Image.fnumber      = 53; %51:55;
XRDIMAGE.Image.numframe     = 1;
XRDIMAGE.Image.numdigs      = 6;
XRDIMAGE.Image.fext         = 'tif';
XRDIMAGE.Image.corrected    = 0;
XRDIMAGE.Image.IsHydra      = 0;    % 0 = Single panel; 1 = GE1; 2 = GE2; 3 = GE3; 4 = GE4;

%%% DARK FILES ONLY USED IF THE IMAGES ARE UNCORRECTED
XRDIMAGE.DarkField.pname    = '/home/beams/S1IDUSER/mnt/s1b/__eval/projects_parkjs/dexela_workflow/shastri_jul19/dunand_AlCe_alloy/data/CeO2_0pt5s';
XRDIMAGE.DarkField.fbase    = 'dark_0pt5s';
XRDIMAGE.DarkField.fnumber  = 61:65;
XRDIMAGE.DarkField.numframe = 1;
XRDIMAGE.DarkField.numdigs  = 6;
XRDIMAGE.DarkField.fext     = 'tif';

%%% CALIBRATION FILE
XRDIMAGE.Calib.pname        = './';
XRDIMAGE.Calib.fbase        = 'CeO2_0pt5s_';
XRDIMAGE.Calib.fnumber      = 53;
XRDIMAGE.Calib.fext         = 'tif';

%%% INSTRUMENT PARAMETERS
XRDIMAGE.Instr.energy       = 71.676;       % keV
XRDIMAGE.Instr.wavelength   = keV2Angstrom(XRDIMAGE.Instr.energy);  % wavelength (Angstrom)
XRDIMAGE.Instr.distance     = 651.224248;  % mm
XRDIMAGE.Instr.centers      = [ 245.854476 , -93.040774 ]; % center offsets x & y (um)
XRDIMAGE.Instr.gammaX       = -0.001182;    % rad
XRDIMAGE.Instr.gammaY       = -0.012298;    % rad
XRDIMAGE.Instr.imrotation   = 0;            %%% NEED TO BE CLARIFIED

% RADIAL CORRECTION
% 0 : no correction
% 1 : constant radial offset
% 2 : PROPOSED BY ISSN 0909-0495 LEE
XRDIMAGE.Instr.dettype  = '2a';

% 0 : []
% 1 : constant value
% 2 : [a1 a2 n1 n2 rhod]
XRDIMAGE.Instr.detpars  = [ ...
    -0.000003 ...
    -0.083582 ...
    -2.496740 ...
    1.171598...
    5687.38873533...
    1.595635 ...
    ];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATERIAL PARAMETERS - CeO2
XRDIMAGE.Material.num       = 1;
XRDIMAGE.Material.lattparms = 5.411651;        % CeO2
XRDIMAGE.Material.structure = 'fcc';
XRDIMAGE.Material.hkls      = load([XRDIMAGE.Material.structure, '.hkls']);

%%% CALCULATE THEORETICAL TTH
[d, th] = PlaneSpacings(XRDIMAGE.Material.lattparms, ...
    'cubic', XRDIMAGE.Material.hkls', ...
    XRDIMAGE.Instr.wavelength);
tth     = 2*th;
d_spacing_range = 0.025;
d_spacing_UB    = (1 + d_spacing_range)*d;
d_spacing_LB    = (1 - d_spacing_range)*d;

tth_UB  = 2.*asind(XRDIMAGE.Instr.wavelength/2)./d_spacing_LB;
tth_LB  = 2.*asind(XRDIMAGE.Instr.wavelength/2)./d_spacing_UB;

XRDIMAGE.Material.tth       = tth;
XRDIMAGE.Material.d_spacing = d;
XRDIMAGE.Material.numpk     = 6;
XRDIMAGE.Material.numbounds = 6;
XRDIMAGE.Material.pkidx     = {...
    [1] [2] [3] [4] [5] [6]
    };
for i = 1:1:XRDIMAGE.Material.numbounds
    XRDIMAGE.Material.pkrange(:,i)  = [ ...
        min(tth_LB(XRDIMAGE.Material.pkidx{i})); ...
        max(tth_UB(XRDIMAGE.Material.pkidx{i})); ...
        ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DATA REDUCTION FLAGS
Analysis_Options.make_polimg    = 1;
Analysis_Options.save_polimg    = 1;
Analysis_Options.fits_spectra   = 1;
Analysis_Options.save_fits      = 1;
Analysis_Options.find_instrpars = 1;
Analysis_Options.save_instrpars = 1;
Analysis_Options.find_detpars	= 1;
Analysis_Options.generateESG    = 0;

%%% PK FITTING OPTIONS
Analysis_Options.PkFuncOptions.pfunc_type	= 'pseudoVoigt';
Analysis_Options.PkFuncOptions.pbkg_order	= 2;
Analysis_Options.PkFitOptimizationOptions   = optimset(...
    'MaxIter', 5e5, ...
    'MaxFunEvals', 3e5);

%%% INSTR OPTIMIZATION OPTIONS
Analysis_Options.InstrPrmFitOptions = optimset(...
    'DerivativeCheck', 'off', ...
    'MaxIter', 1e5, ...
    'MaxFunEvals', 3e5, ...
    'TypicalX',[100 -100 682 0.1 0.1 XRDIMAGE.Instr.detpars], ...
    'Display','final');

fname_pattern   = sprintf('%%s_%%0%dd.%%s', XRDIMAGE.Image.numdigs);
    
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
        pfname  = fullfile(XRDIMAGE.Image.pname, fname);
        disp('###########################')
        disp(sprintf('Looking at %s', pfname))
        disp('###########################')
        
        pfname_polimage = [pfname, '.polimg.mat'];
        
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
        
        if Analysis_Options.save_polimg
            disp(sprintf('Saving polimg for %s', pfname))
            save(pfname_polimage, 'polimg', 'XRDIMAGE')
        else
            disp(sprintf('Not saving polimg for %s', pfname))
        end
        
        figure(2)
        subplot(1,2,1)
        imagesc(log(polimg.intensity)), axis square tight
        hold off
        
        subplot(1,2,2)
        plot(polimg.radius, polimg.intensity)
        hold off
        disp(' ')
    end
end

if Analysis_Options.fits_spectra
    for i = 1:1:numimg    
        fname	= sprintf(fname_pattern, XRDIMAGE.Image.fbase, XRDIMAGE.Image.fnumber, XRDIMAGE.Image.fext);
        pfname  = fullfile(XRDIMAGE.Image.pname, fname);
        
        pfname_polimage = [pfname, '.polimg.mat'];
        pfname_pkfit    = [pfname, '.pkfit.mat'];
        
        polimg  = load(pfname_polimage);
        polimg  = polimg.polimg;
        
        figure(2)
        subplot(1,2,1)
        imagesc(log(polimg.intensity)), axis square tight
        hold off
        
        subplot(1,2,2)
        plot(polimg.radius, polimg.intensity)
        hold off
        
        disp('###########################')
        disp(sprintf('Fitting peaks in %s', pfname))
        disp('###########################')
        for j = 1:1:XRDIMAGE.CakePrms.bins(1)
            disp(sprintf('Looking at azimuthal bin %d of %d\n', j, XRDIMAGE.CakePrms.bins(1)))
            
            x   = polimg.radius(j,:);
            x   = Pixel2mm(x, XRDIMAGE.Instr.pixelsizeHorz);  % CONVERT TO MM FROM PIXELS
            y   = polimg.intensity(j,:);
            
            figure(11)
            subplot(1,2,1)
            plot(x, y, 'k.')
            hold on
            plot(XRDIMAGE.Instr.distance*tand(tth), mean(y), 'g^')
            axis([min(x) max(x) 0 max(y)+100])
            xlabel('radial distance (mm)')
            ylabel('intensity (arb. units)')
            title(['bin number : ', num2str(j)])
            
            for k = 1:1:XRDIMAGE.Material.numpk
                disp(sprintf('Looking at peak number %d of %d', k, XRDIMAGE.Material.numpk))
                if j == 1
                    pkrange = XRDIMAGE.Material.pkrange(:,k);
                    pkrange = XRDIMAGE.Instr.distance.*tand(pkrange);
                    
                    idx = find(x >= pkrange(1) & x <= pkrange(2));
                    xr  = x(idx)';
                    yr  = y(idx)';
                    
                    pr0 = [...
                        max(yr)/5 ...
                        0.5 ...
                        0.15 ...
                        XRDIMAGE.Instr.distance*tand(tth(XRDIMAGE.Material.pkidx{k})) ...
                        0 ...
                        0];
                else
                    pkrange = [pkfit.rho(j-1,k)-2.5 pkfit.rho(j-1,k)+2.5];
                    idx = find(x >= pkrange(1) & x <= pkrange(2));
                    xr  = x(idx)';
                    yr  = y(idx)';
                    
                    pr0 = pr(k,:);
                end
                
                pkpars.pfunc_type   = Analysis_Options.PkFuncOptions.pfunc_type;
                pkpars.pbkg_order   = Analysis_Options.PkFuncOptions.pbkg_order;
                pkpars.xdata        = xr;
                
                [pr(k,:), rsn, ~, ef]   = lsqcurvefit(@pfunc_switch, pr0, pkpars, yr, ...
                    [], [], Analysis_Options.PkFitOptimizationOptions);
                
                y0	= pfunc_switch(pr0, pkpars);
                yf	= pfunc_switch(pr(k,:), pkpars);
                
                figure(11)
                subplot(1,2,1)
                plot(xr, yr, 'b.')
                plot(xr, y0, 'r:')
                plot(xr, yf, 'g-')
                
                subplot(1,2,2)
                plot(xr, yr, 'b.')
                hold on
                plot(xr, y0, 'r:')
                plot(xr, yf, 'g-')
                xlabel('radial distance (mm)')
                ylabel('intensity (arb. units)')
                title(['peak number : ', num2str(k)])
                hold off
                
                %%% MAPPING NEEDS TO BE UPDATED WITH A SWITCH
                pro  = pkfit_MapResult(pkpars, pr(k,:));
                
                pkfit.amp(j,k)  = pro(1);
                pkfit.fwhm{j,k} = pro(2:3);
                pkfit.mix{j,k}  = pro(4:5);
                pkfit.rho(j,k)  = pro(6);
                pkfit.bkg{j,k}  = pro(7:end);
                pkfit.rsn(j,k)  = rsn;
                pkfit.ef(j,k)   = ef;
                pkfit.rwp(j,k)  = ErrorRwp(yr, yf);
            end
            figure(11)
            subplot(1,2,1)
            hold off
        end
        
        if Analysis_Options.save_fits
            disp('###########################')
            disp(sprintf('Saving peak fits in %s\n', pfname_pkfit))
            save(pfname_pkfit, 'pkfit')
        else
            disp('###########################')
            disp(sprintf('Not saving peak fits for %s\n', pfname{i,1}))
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% APPLY/FIND GEOMETRICAL MODEL
if Analysis_Options.find_instrpars
    for i = 1:1:numimg
        fname	= sprintf(fname_pattern, XRDIMAGE.Image.fbase, XRDIMAGE.Image.fnumber, XRDIMAGE.Image.fext);
        pfname  = fullfile(XRDIMAGE.Image.pname, fname);
        
        pfname_polimage = [pfname, '.polimg.mat'];
        pfname_pkfit    = [pfname, '.pkfit.mat'];
        pfname_instr    = [pfname, '.instr.mat'];
        
        disp('###########################')
        disp(sprintf('Looking at %s to find instrument parameters.', pfname))
        disp('###########################')
        
        disp(sprintf('Loading polimg in %s\n', pfname_polimage))
        polimg  = load(pfname_polimage);
        polimg  = polimg.polimg;
        
        disp(sprintf('Loading peak fits in %s\n', pfname_pkfit))
        load(pfname_pkfit)
        
        p0  = [...
            XRDIMAGE.Instr.centers, ...
            XRDIMAGE.Instr.distance, ...
            XRDIMAGE.Instr.gammaX, ...
            XRDIMAGE.Instr.gammaY, ...
            XRDIMAGE.Instr.detpars];
        
        GeomModelParams.pkidx           = [XRDIMAGE.Material.pkidx{:}]';
        GeomModelParams.tth             = XRDIMAGE.Material.tth;
        GeomModelParams.azim            = XRDIMAGE.CakePrms.azim;
        GeomModelParams.rho             = pkfit.rho';
        GeomModelParams.dettype         = XRDIMAGE.Instr.dettype;
        GeomModelParams.DistortParams0  = XRDIMAGE.Instr.detpars;
        GeomModelParams.find_detpars    = Analysis_Options.find_detpars;
        
        Instr	= XRDIMAGE.Instr;
        
        dtth0   = ApplyGeometricModel(p0, GeomModelParams);
        tth0    = tth([XRDIMAGE.Material.pkidx{:}])';
        tth0    = repmat(tth0, 1, size(dtth0, 2));
        strain0 = sind(tth0)./sind(tth0 - dtth0) - 1;
        
        ydata   = zeros(XRDIMAGE.Material.numpk, XRDIMAGE.CakePrms.bins(1));
        p       = lsqcurvefit(@ApplyGeometricModel, p0, GeomModelParams, ydata, [], [], Analysis_Options.InstrPrmFitOptions);
        
        dtth    = ApplyGeometricModel(p, GeomModelParams);
        strain  = sind(tth0)./sind(tth0 - dtth) - 1;
        
        %%% ASSIGN NEW INSTRUMENT PARAMETERS USING OPTIMIZATION RESULTS
        Instr.centers   = p(1:2);
        Instr.distance  = p(3);
        Instr.gammaX    = p(4);
        Instr.gammaY    = p(5);
        Instr.detpars   = p(6:end);
        
        disp('Instrument parameter optimization results')
        disp('Update parameters in XRDIMAGE.Instr variable accordingly')
        disp(sprintf('Instr.centers  : %f , %f', p(1), p(2)))
        disp(sprintf('Instr.distance : %f', p(3)))
        disp(sprintf('Instr.gammaX   : %f', p(4)))
        disp(sprintf('Instr.gammaY   : %f', p(5)))
        disp(sprintf('Detector distortion prm : %f\n', p(6:end)))
        
        figure(100)
        subplot(1,2,1)
        imagesc(strain0')
        colorbar vert
        title('pseudo-strain due to p0')
        xlabel('hkl id')
        ylabel('azimuthal bin number')
        
        subplot(1,2,2)
        imagesc(strain')
        colorbar vert
        title('pseudo-strain due to p')
        xlabel('hkl id')
        ylabel('azimuthal bin number')
        
        mapped_tth  = GeometricModelXRDSwitch(Instr, polimg);
        polimg.mapped_tth_for_intensity = mapped_tth;
        
        [tth_grid, d_grid, intensity_in_tth_grid]   = MapIntensityToTThGrid(XRDIMAGE, polimg);
        polimg.tth_grid                 = tth_grid;
        polimg.d_grid                   = d_grid;
        polimg.intensity_in_tth_grid    = intensity_in_tth_grid;
        
        figure(101)
        imagesc(log(abs(polimg.intensity_in_tth_grid))), axis square tight
        title('Caked image // radial position is corrected')
        hold off
        
        disp('###########################')
        disp(sprintf('mean pseudo-strain using p0 : %f', mean(abs(strain0(:)))))
        disp(sprintf('mean pseudo-strain using p  : %f\n', mean(abs(strain(:)))))
        disp(sprintf('std pseudo-strain using p0 : %f', std(strain0(:))))
        disp(sprintf('std pseudo-strain using p  : %f\n', std(strain(:))))
        
        if Analysis_Options.save_instrpars
            disp('###########################')
            disp(sprintf('Saving optimized innstrument parameters in %s\n', pfname_instr))
            save(pfname_instr, 'Instr')
        else
            disp('###########################')
            disp(sprintf('NOT saving optimized instrument parameters for %s\n', pfname{i,1}))
        end
    end
end
