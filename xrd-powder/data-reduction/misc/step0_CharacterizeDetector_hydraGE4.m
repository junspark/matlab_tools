clear all
close all
clc

%%% INPUT PARAMETERS
XRDIMAGE.Image.pname        = 'S:\park_jul16\ge2\';
XRDIMAGE.Image.fbase        = 'CeO2_3s_';
XRDIMAGE.Image.fnumber      = 955;
XRDIMAGE.Image.numframe     = 5;
XRDIMAGE.Image.numdigs      = 5;
XRDIMAGE.Image.fext         = 'ge4';
XRDIMAGE.Image.corrected    = 0;
XRDIMAGE.Image.IsHydra      = 1;    % 0 = Single panel; 1 = GE1; 2 = GE2; 3 = GE3; 4 = GE4;

%%% DARK FILES ONLY USED IF THE IMAGES ARE UNCORRECTED
XRDIMAGE.DarkField.pname    = 'S:\park_jul16\ge1\';
XRDIMAGE.DarkField.fbase    = 'dark_3s_';
XRDIMAGE.DarkField.fnumber  = 732;
XRDIMAGE.DarkField.numframe = 1;
XRDIMAGE.DarkField.numdigs  = 5;
XRDIMAGE.DarkField.fext     = 'ge4';

XRDIMAGE.Calib.pname        = 'S:\park_jul16\ge2\';
XRDIMAGE.Calib.fbase        = 'CeO2_3s_';
XRDIMAGE.Calib.fnumber      = 734;
XRDIMAGE.Calib.fext         = 'ge4';

%%% INSTRUMENT PARAMETERS
XRDIMAGE.Instr.energy       = 71.676;       % keV
XRDIMAGE.Instr.wavelength   = keV2Angstrom(XRDIMAGE.Instr.energy);  % wavelength (Angstrom)
XRDIMAGE.Instr.distance     = 1.89516943e+03;     % mm
XRDIMAGE.Instr.centers      = [ 950.294775 , 720.273581 ]; % center offsets x & y (um)
XRDIMAGE.Instr.gammaX       = -0.014297;    % rad
XRDIMAGE.Instr.gammaY       = 0.010257;    % rad
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
    0.000440 ...
    0.000104 ...
    -0.454242 ...
    -1.516900 ...
    3099.498067 ...
    -3.199170 ...
    ];

%%% CAKE PARAMETERS
XRDIMAGE.CakePrms.bins(1)   = 10;               % number of azimuthal bins over angular range defined by XRDIMAGE.CakePrms.sector(1) and XRDIMAGE.CakePrms.sector(2)
XRDIMAGE.CakePrms.bins(2)   = 3000;             % number of radial bins over radial range defined by XRDIMAGE.CakePrms.sector(3) and XRDIMAGE.CakePrms.sector(4)
XRDIMAGE.CakePrms.bins(3)   = 5;               % number of angular bins

XRDIMAGE.CakePrms.origin(1) = 4.90258858e+02;         % apparent X center in hexrd (mm)
XRDIMAGE.CakePrms.origin(2) = -1.52089437e+01;        % apparent Y center in hexrd (mm)
XRDIMAGE.CakePrms.origin(1) = XRDIMAGE.CakePrms.origin(1)/XRDIMAGE.Instr.pixelsizeHorz;         % convert to pixels
XRDIMAGE.CakePrms.origin(2) = XRDIMAGE.CakePrms.origin(2)/XRDIMAGE.Instr.pixelsizeVert;         % convert to pixels
XRDIMAGE.CakePrms.origin(2) = 2048-XRDIMAGE.CakePrms.origin(2); %%% CONVERT TO IMAGE COORDINATES

XRDIMAGE.CakePrms.sector(1) = 200;      % start azimuth (min edge of bin) in degrees
XRDIMAGE.CakePrms.sector(2) = 240;      % stop  azimuth (max edge of bin) in degrees
XRDIMAGE.CakePrms.sector(3) = 600;      % start radius (min edge of bin) in pixels
XRDIMAGE.CakePrms.sector(4) = 2300;     % stop  radius (max edge of bin) in pixels

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
d_spacing_range = 0.05;
d_spacing_UB    = (1 + d_spacing_range)*d;
d_spacing_LB    = (1 - d_spacing_range)*d;

tth_UB  = 2.*asind(XRDIMAGE.Instr.wavelength/2)./d_spacing_LB;
tth_LB  = 2.*asind(XRDIMAGE.Instr.wavelength/2)./d_spacing_UB;

XRDIMAGE.Material.tth       = tth;
XRDIMAGE.Material.d_spacing = d;

% XRDIMAGE.Material.numpk     = 10;
% XRDIMAGE.Material.numbounds = 10;
% XRDIMAGE.Material.pkidx     = {...
%     [3] [4] [5] [6] [7] [8] [9] [12] [13] [16]
%     };

XRDIMAGE.Material.numpk     = 5;
XRDIMAGE.Material.numbounds = 5;
XRDIMAGE.Material.pkidx     = {...
    [3] [6] [9] [12] [16]
    };

% XRDIMAGE.Material.numpk     = 1;
% XRDIMAGE.Material.numbounds = 1;
% XRDIMAGE.Material.pkidx     = {...
%     [1]
%     };

% XRDIMAGE.Material.numpk     = 2;
% XRDIMAGE.Material.numbounds = 2;
% XRDIMAGE.Material.pkidx     = {...
%     [1] [2]
%     };

% XRDIMAGE.Material.numpk     = 3;
% XRDIMAGE.Material.numbounds = 3;
% XRDIMAGE.Material.pkidx     = {...
%     [7] [9] [10]
%     };

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

Analysis_Options.InstrPrmFitOptions = optimset(...
    'DerivativeCheck', 'off', ...
    'MaxIter', 1e5, ...
    'MaxFunEvals', 3e5, ...
    'TypicalX',[100 -100 100 0.1 0.1 XRDIMAGE.Instr.detpars], ...
    'Display','iter');

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
        
        figure(1)
        hold off
        imagesc(rot90(imgi,1))
        caxis([-10 3000])
        axis equal
        colorbar vert
        hold on
        xlabel('X_L (pixels)')
        ylabel('Y_L (pixels)')
        title('Ensure that image matches the coordinate system')
        
        %%% POLAR REBINNING
        polimg  = PolarBinXRD(DetectorMesh, ...
            XRDIMAGE.Instr, ...
            XRDIMAGE.CakePrms, ...
            imgi);
        
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
                        2e2];
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
                
                %%% NEEDS TO BE ADAPTIVE FOR PEAK FUNCTION TYPE
                pLB = [0 0 0 -inf -inf -inf];
                pUB = [inf inf 1 inf inf inf];
                
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
        
        XRDIMAGE.Instr  = Instr;
        
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
