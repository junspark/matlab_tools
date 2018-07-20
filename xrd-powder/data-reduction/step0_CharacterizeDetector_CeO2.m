%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instrument parameter optimization results
% Update parameters in XRDIMAGE.Instr variable accordingly
% Instr.centers  : -146.440630 , -366.383718
% Instr.distance : 1001.306067
% Instr.gammaX   : -0.002070
% Instr.gammaY   : -0.001311
% Detector distortion prm : 0.000176
% Detector distortion prm : -0.001961
% Detector distortion prm : 1.929915
% Detector distortion prm : 2.282699
% Detector distortion prm : 417.585333
% Detector distortion prm : 2.122985
% 
% GeometricModelXRD2a
% ###########################
% mean pseudo-strain using p0 : 0.000246
% mean pseudo-strain using p  : 0.000022
% 
% std pseudo-strain using p0 : 0.000293
% std pseudo-strain using p  : 0.000030
% ###########################
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%%% INPUT PARAMETERS
XRDIMAGE.Image.pname        = '/home/beams/S1IDUSER/mnt/s1b/__eval/projects_parkjs/dexela_calibration/clausen_jul18/ge5_in_e';
XRDIMAGE.Image.fbase        = 'CeO2_0pt3s_';
XRDIMAGE.Image.fnumber      = 53; % 4116 / 4117
XRDIMAGE.Image.numframe     = 5;
XRDIMAGE.Image.numdigs      = 6;
XRDIMAGE.Image.fext         = 'ge5';
XRDIMAGE.Image.corrected    = 0;
XRDIMAGE.Image.IsHydra      = 0;    % 0 = Single panel; 1 = GE1; 2 = GE2; 3 = GE3; 4 = GE4;

%%% DARK FILES ONLY USED IF THE IMAGES ARE UNCORRECTED
XRDIMAGE.DarkField.pname    = '/home/beams/S1IDUSER/mnt/s1b/__eval/projects_parkjs/dexela_calibration/clausen_jul18/ge5_in_e';
XRDIMAGE.DarkField.fbase    = 'Dark_0pt3s_';
XRDIMAGE.DarkField.fnumber  = 54;
XRDIMAGE.DarkField.numframe = 5;
XRDIMAGE.DarkField.numdigs  = 6;
XRDIMAGE.DarkField.fext     = 'ge5';

XRDIMAGE.Calib.pname        = '/home/beams/S1IDUSER/mnt/s1b/__eval/projects_parkjs/dexela_calibration/clausen_jul18/ge5_in_e';
XRDIMAGE.Calib.fbase        = 'CeO2_0pt3s_';
XRDIMAGE.Calib.fnumber      = 53;
XRDIMAGE.Calib.fext         = 'ge5';

%%% INSTRUMENT PARAMETERS
XRDIMAGE.Instr.energy       = 71.676;       % keV
XRDIMAGE.Instr.wavelength   = keV2Angstrom(XRDIMAGE.Instr.energy);  % wavelength (Angstrom)
XRDIMAGE.Instr.distance     = 1390.815364;  % mm
XRDIMAGE.Instr.centers      = [ -285.632569 , 2.094275 ]; % center offsets x & y (um)
XRDIMAGE.Instr.gammaX       = -0.003136;    % rad
XRDIMAGE.Instr.gammaY       = -0.008367;    % rad
XRDIMAGE.Instr.detsizeHorz  = 409.6;    % mm
XRDIMAGE.Instr.detsizeVert  = 409.6;    % mm
XRDIMAGE.Instr.pixelsizeHorz    = 0.2;          % mm
XRDIMAGE.Instr.pixelsizeVert    = 0.2;          % mm
XRDIMAGE.Instr.numpixelsHorz    = XRDIMAGE.Instr.detsizeHorz/XRDIMAGE.Instr.pixelsizeHorz;   % total number of rows in the full image
XRDIMAGE.Instr.numpixelsVert    = XRDIMAGE.Instr.detsizeVert/XRDIMAGE.Instr.pixelsizeVert;   % total number of rows in the full image
XRDIMAGE.Instr.imrotation   = 0;

% RADIAL CORRECTION
% 0 : no correction
% 1 : constant radial offset
% 2 : PROPOSED BY ISSN 0909-0495 LEE
XRDIMAGE.Instr.dettype  = '2a';

% 0 : []
% 1 : constant value
% 2 : [a1 a2 n1 n2 rhod]
XRDIMAGE.Instr.detpars  = [ ...
    -0.00001756 ...
    -0.00002753 ...
    0.04495597 ...
    0.01995934...
    5.38873533...
    0.01873452 ...
    ]*1e2;

%%% CAKE PARAMETERS
XRDIMAGE.CakePrms.bins(1)   = 72;           % number of azimuthal bins
XRDIMAGE.CakePrms.bins(2)   = 1000;         % number of radial bins
XRDIMAGE.CakePrms.bins(3)   = 70;            % number of angular bins

XRDIMAGE.CakePrms.origin(1) = 1024.190;     % apparent X center in pixels // THIS IS WHAT YOU SEE ON FIGURE 1
XRDIMAGE.CakePrms.origin(2) = 1025.110;     % apparent Y center in pixels // THIS IS WHAT YOU SEE ON FIGURE 1
XRDIMAGE.CakePrms.origin(2) = XRDIMAGE.Instr.numpixelsVert-XRDIMAGE.CakePrms.origin(2); %%% CONVERT TO IMAGE COORDIANTES

XRDIMAGE.CakePrms.sector(1) = -360/XRDIMAGE.CakePrms.bins(1)/2;     % start azimuth (min edge of bin) in degrees
XRDIMAGE.CakePrms.sector(2) = 360-360/XRDIMAGE.CakePrms.bins(1)/2;  % stop  azimuth (max edge of bin) in degrees
XRDIMAGE.CakePrms.sector(3) = 150;  % start radius (min edge of bin) in pixels
XRDIMAGE.CakePrms.sector(4) = 1000;  % stop  radius (max edge of bin) in pixels

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
d_spacing_range = 0.01;
d_spacing_UB    = (1 + d_spacing_range)*d;
d_spacing_LB    = (1 - d_spacing_range)*d;

tth_UB  = 2.*asind(XRDIMAGE.Instr.wavelength/2)./d_spacing_LB;
tth_LB  = 2.*asind(XRDIMAGE.Instr.wavelength/2)./d_spacing_UB;

XRDIMAGE.Material.tth       = tth;
XRDIMAGE.Material.d_spacing = d;
XRDIMAGE.Material.numpk     = 4;
XRDIMAGE.Material.numbounds = 4;
XRDIMAGE.Material.pkidx     = {...
    [2] [3] [4] [5]
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
    'MaxFunEvals',3e5);

%%% INSTR OPTIMIZATION OPTIONS
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
    pfname  = GenerateGEpfname(XRDIMAGE.DarkField);
    bg      = NreadGE(pfname{1,1}, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GENERATE MESH FOR INTEGRATION
%%% IF POLIMG NEEDS TO BE GENERATED
if Analysis_Options.make_polimg
    if ~XRDIMAGE.CakePrms.fastint
        DetectorMesh    = BuildMeshDetector(XRDIMAGE.Instr.numpixelsHorz, XRDIMAGE.Instr.numpixelsVert, XRDIMAGE.CakePrms);
    else
        DetectorMesh    = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOAD XRD IMAGES & GENERATE POLIMG IF NECESSARY
pfname  = GenerateGEpfname(XRDIMAGE.Image);
numimg  = length(pfname);
if Analysis_Options.make_polimg
    for i = 1:1:numimg
        disp('###########################')
        disp(sprintf('Looking at %s', pfname{i,1}))
        disp('###########################')
        
        pfname_polimage = [pfname{i,1}, '.polimg.mat'];
        
        if XRDIMAGE.Image.corrected
            imgi    = ReadSUM(pfname{i,1});
        else
            imgi    = bg.*0;
            for j = 1:1:XRDIMAGE.Image.numframe
                imgj    = NreadGE(pfname{i,1}, j);
                imgi    = imgi + imgj;
            end
            imgi    = imgi - bg.*XRDIMAGE.Image.numframe;
        end
        
        %%% SYNTHETIC IMAGE
%         xxx = (1:1:2048)./204.8;
%         yyy = (1:1:2048)./(204.8*2);
%         imgi    = xxx'*yyy;
%         imgi(1000:1200,:)   = 27;
        
        figure(1)
        hold off
        imagesc(rot90(imgi,1)) %% DO NOT CHANGE
        % caxis([-10 3000])
        axis equal tight
        colorbar vert
        hold on
        xlabel('X_L (pixels)')
        ylabel('Y_L (pixels)')
        title('Ensure that image matches the coordinate system')
        text(0, 0, 'TO')
        text(2048, 0, 'TI')
        text(0, 0, 'TO')
        text(0, 2048, 'BO')
        text(2048, 2048, 'BI')
        
        %%% POLAR REBINNING
        polimg	= PolarBinXRD(DetectorMesh, ...
            XRDIMAGE.Instr, ...
            XRDIMAGE.CakePrms, ...
            imgi, 'PlotProgress', 'off');
        
        if Analysis_Options.save_polimg
            disp(sprintf('Saving polimg for %s', pfname{i,1}))
            save(pfname_polimage, 'polimg', 'XRDIMAGE')
        else
            disp(sprintf('Not saving polimg for %s', pfname{i,1}))
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
        pfname_polimage = [pfname{i,1}, '.polimg.mat'];
        pfname_pkfit    = [pfname{i,1}, '.pkfit.mat'];
        
        load(pfname_polimage)
        figure(2)
        subplot(1,2,1)
        imagesc(log(polimg.intensity)), axis square tight
        hold off
        
        subplot(1,2,2)
        plot(polimg.radius, polimg.intensity)
        hold off
        
        disp('###########################')
        disp(sprintf('Fitting peaks in %s', pfname{i,1}))
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
        pfname_polimage = [pfname{i,1}, '.polimg.mat'];
        pfname_pkfit    = [pfname{i,1}, '.pkfit.mat'];
        pfname_instr    = [pfname{i,1}, '.instr.mat'];
        
        disp('###########################')
        disp(sprintf('Looking at %s to find instrument parameters.', pfname{i,1}))
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
        
        Instr           = XRDIMAGE.Instr;
        
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
