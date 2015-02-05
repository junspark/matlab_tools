clear all
close all
clc

%%% INPUT PARAMETERS
XRDIMAGE.DarkField.pname    = 'O:\park_nov14';
XRDIMAGE.DarkField.fbase    = 'dark_1.5s_';
XRDIMAGE.DarkField.fnumber  = 338;
XRDIMAGE.DarkField.numframe = 1;
XRDIMAGE.DarkField.numdigs  = 5;
XRDIMAGE.DarkField.fext     = 'ge3';

XRDIMAGE.Image.pname        = 'O:\park_nov14';
XRDIMAGE.Image.fbase        = 'CeO2_';
XRDIMAGE.Image.fnumber      = 746; %
XRDIMAGE.Image.numframe     = 1;
XRDIMAGE.Image.numdigs      = 5;
XRDIMAGE.Image.fext         = 'ge3.sum';
XRDIMAGE.Image.corrected    = 1;

XRDIMAGE.Calib.pname        = 'O:\park_nov14';
XRDIMAGE.Calib.fbase        = 'CeO2_1.5s_';
XRDIMAGE.Calib.fnumber      = 746;

%%% INSTRUMENT PARAMETERS
XRDIMAGE.Instr.energy       = 88.005;       % keV
XRDIMAGE.Instr.wavelength   = keV2Angstrom(XRDIMAGE.Instr.energy);  % wavelength (Angstrom)
XRDIMAGE.Instr.detectorsize = 409.6;        % mm
XRDIMAGE.Instr.pixelsize    = 0.2;          % mm
XRDIMAGE.Instr.distance     = 1650.695095;   % mm
XRDIMAGE.Instr.centers      = [611.898998 , 468.762851]; % center offsets x & y (um)
XRDIMAGE.Instr.gammaX       = 0.003947;    % rad
XRDIMAGE.Instr.gammaY       = -0.000077;    % rad
XRDIMAGE.Instr.numpixels    = XRDIMAGE.Instr.detectorsize/XRDIMAGE.Instr.pixelsize;   % total number of rows in the full image

% RADIAL CORRECTION
% 0 : no correction
% 1 : constant radial offset
% 2 : PROPOSED BY ISSN 0909-0495 LEE
XRDIMAGE.Instr.dettype  = '2a';

% 0 : []
% 1 : constant value
% 2 : [a1 a2 n1 n2 rhod]
XRDIMAGE.Instr.detpars  = [ ...
    -0.000018  ...
    -0.000259  ...
    0.211526  ...
    3.214555  ...
    210.611340  ...
    2.306058 ...
    ];

%%% CAKE PARAMETERS
XRDIMAGE.CakePrms.bins(1)   = 72;           % number of azimuthal bins
XRDIMAGE.CakePrms.bins(2)   = 2000;         % number of radial bins
XRDIMAGE.CakePrms.bins(3)   = 20;            % number of angular bins
XRDIMAGE.CakePrms.origin(1) = 1038.159;         % x center in pixels, fit2d Y 
XRDIMAGE.CakePrms.origin(2) = 1013.426;         % y center in pixels, fit2d X
XRDIMAGE.CakePrms.sector(1) = -360/XRDIMAGE.CakePrms.bins(1)/2;     % start azimuth (min edge of bin) in degrees
XRDIMAGE.CakePrms.sector(2) = 360-360/XRDIMAGE.CakePrms.bins(1)/2;  % stop  azimuth (max edge of bin) in degrees
XRDIMAGE.CakePrms.sector(3) = 200;  % start radius (min edge of bin) in pixels
XRDIMAGE.CakePrms.sector(4) = 950;  % stop  radius (max edge of bin) in pixels
XRDIMAGE.CakePrms.azim      = 0:360/XRDIMAGE.CakePrms.bins(1):XRDIMAGE.CakePrms.sector(2);
XRDIMAGE.CakePrms.fastint   = 1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% MATERIAL PARAMETERS - CeO2
XRDIMAGE.Material.num       = 1;
XRDIMAGE.Material.lattparms = 5.411651;        % CeO2
XRDIMAGE.Material.structure = 'fcc';
XRDIMAGE.Material.numpk     = 7;
XRDIMAGE.Material.pkidx     = {...
    [1] [2] [3] [4] [5] [6] [7]
    };

XRDIMAGE.Material.pkrange    = [...
    2.5337    2.9335    4.1698    4.8986    5.1188   5.9391    6.4552    6.6544; ...
    2.6637    3.0635    4.2898    5.0186    5.2488   6.0291    6.5752    6.6944; ...
    ];

XRDIMAGE.Material.pkbck     = 2;
XRDIMAGE.Material.pkfunc    = 4;
XRDIMAGE.Material.hkls      = load([XRDIMAGE.Material.structure, '.hkls']);

%%% CALCULATE THEORETICAL TTH
[d, th] = PlaneSpacings(XRDIMAGE.Material.lattparms, ...
    'cubic', XRDIMAGE.Material.hkls', ...
    XRDIMAGE.Instr.wavelength);
tth     = 2*th;

XRDIMAGE.Material.tth       = tth;
XRDIMAGE.Material.d_spacing = d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% MATERIAL PARAMETERS - LaB6
% XRDIMAGE.Material.num       = 1;
% XRDIMAGE.Material.lattparms = 4.1569162;        % LaB6
% XRDIMAGE.Material.structure = 'simplecubic';
% XRDIMAGE.Material.numpk     = 8;
% XRDIMAGE.Material.pkrange    = [...
%     2.7  3.8 4.7 6.1 6.7 8.2 8.7 9.1; ...
%     2.95 4.1 5.0 6.7 7.0 8.5 9.0 9.4; ...
%     ];
% XRDIMAGE.Material.pkidx     = {...
%     [1] [2] [3] [5] [6] [8] [9] [10]
%     };
% XRDIMAGE.Material.pkbck     = 2;
% XRDIMAGE.Material.pkfunc    = 4;
% XRDIMAGE.Material.hkls      = load([XRDIMAGE.Material.structure, '.hkls']);
% 
% %%% CALCULATE THEORETICAL TTH
% [d, th] = PlaneSpacings(XRDIMAGE.Material.lattparms, ...
%     'cubic', XRDIMAGE.Material.hkls', ...
%     XRDIMAGE.Instr.wavelength);
% tth     = 2*th;
% 
% XRDIMAGE.Material.tth       = tth;
% XRDIMAGE.Material.d_spacing = d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DATA REDUCTION FLAGS
Analysis_Options.make_polimg    = 1;
Analysis_Options.save_polimg    = 1;
Analysis_Options.fits_spectra   = 1;
Analysis_Options.save_fits      = 1;
Analysis_Options.find_instrpars = 1;
Analysis_Options.save_instrpars = 1;
Analysis_Options.find_detpars	= 1;

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
    DetectorMesh    = BuildMeshDetector(XRDIMAGE.Instr.numpixels);
end

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
%%% LOAD XRD IMAGES & GENERATE POLIMG
pfname  = GenerateGEpfname(XRDIMAGE.Image);
numimg  = length(pfname);
if Analysis_Options.make_polimg
    for i = 1:1:numimg
        disp('###########################')
        disp(sprintf('Looking at %s', pfname{i,1}))
        disp('###########################')
        
        pfname_polimage = [pfname{i,1}, '.polimg.mat'];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% POLAR REBINNING IF NECESSARY
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
        
        % NOTE ROTATION IS NEED TO ENSURE PROPER REBINNING
        imgi    = rot90(imgi);
        
        figure(1)
        imagesc(imgi)
        axis equal tight
        colorbar vert
        xlabel('X_{L}')
        ylabel('Y_{L}')
        hold off
        
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
        disp(sprintf('Fitting peaks in %s', pfname{i,1}))
        disp('###########################')
        
        for j = 1:1:XRDIMAGE.CakePrms.bins(1)
            disp(sprintf('Looking at azimuthal bin %d of %d\n', j, XRDIMAGE.CakePrms.bins(1)))
            
            x   = polimg.radius(j,:);
            x   = Pixel2mm(x, XRDIMAGE.Instr.pixelsize);  % CONVERT TO MM FROM PIXELS
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
                        2e3];
                else
                    pkrange = [pkfit.rho(j-1,k)-2.5 pkfit.rho(j-1,k)+2.5];
                    idx = find(x >= pkrange(1) & x <= pkrange(2));
                    xr  = x(idx)';
                    yr  = y(idx)';
                    
                    pr0 = [...
                        pkfit.amp(j-1,k) ...
                        pkfit.mix(j-1,k) ...
                        pkfit.fwhm(j-1,k) ...
                        pkfit.rho(j-1,k) ...
                        pkfit.bkg{j-1,k}];
                end
                
                y0  = pfunc(pr0,xr);
                [pr, rsn, ~, ef]    = lsqcurvefit(@pfunc, pr0, xr, yr, ...
                    [], [], Analysis_Options.PkFitOptions);
                yf  = pfunc(pr,xr);
                
                figure(11)
                subplot(1,2,1)
                plot(xr, yr, 'b.')
                plot(xr, y0, 'r-')
                plot(xr, yf, 'g-')
                
                subplot(1,2,2)
                plot(xr, yr, 'b.')
                hold on
                plot(xr, y0, 'r-')
                plot(xr, yf, 'g-')
                xlabel('radial distance (mm)')
                ylabel('intensity (arb. units)')
                title(['peak number : ', num2str(k)])
                hold off
                
                pkfit.amp(j,k)  = pr(1);
                pkfit.mix(j,k)  = pr(2);
                pkfit.fwhm(j,k) = pr(3);
                pkfit.rho(j,k)  = pr(4);
                pkfit.bkg{j,k}  = pr(5:end);
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
        
        disp(sprintf('Loading peak fits in %s\n', pfname_polimage))
        polimg  = load(pfname_polimage);
        polimg  = polimg.polimg;
        
        disp(sprintf('Loading peak fits in %s\n', pfname_pkfit))
        pkfit   = load(pfname_pkfit);
        pkfit   = pkfit.pkfit;
        
        p0  = [...
            XRDIMAGE.Instr.centers, ...
            XRDIMAGE.Instr.distance, ...
            XRDIMAGE.Instr.gammaX, ...
            XRDIMAGE.Instr.gammaY, ...
            XRDIMAGE.Instr.detpars];
        
        dtth0   = ApplyGeometricModel(p0);
        tth0    = tth([XRDIMAGE.Material.pkidx{:}])';
        tth0    = repmat(tth0, 1, size(dtth0, 2));
        strain0 = sind(tth0)./sind(tth0 - dtth0) - 1;
        
        p	= lsqnonlin(@ApplyGeometricModel, p0, [], [], Analysis_Options.InstrPrmFitOptions);
        
        disp('Instrument parameter optimization results')
        disp('Update parameters accordingly')
        disp(sprintf('XRDIMAGE.Instr.centers  : %f , %f', p(1), p(2)))
        disp(sprintf('XRDIMAGE.Instr.distance : %f', p(3)))
        disp(sprintf('XRDIMAGE.Instr.gammaX   : %f', p(4)))
        disp(sprintf('XRDIMAGE.Instr.gammaY   : %f', p(5)))
        disp(sprintf('Detector distortion prm : %f\n', p(6:end)))
        
        dtth    = ApplyGeometricModel(p);
        strain  = sind(tth0)./sind(tth0 - dtth) - 1;
        
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
        
        Data    = cell(1, XRDIMAGE.CakePrms.bins(1));
        for ii=1:1:XRDIMAGE.CakePrms.bins(1)
            Data{ii}    = [XRDIMAGE.Instr.pixelsize*polimg.radius(ii,:)' polimg.intensity(ii,:)'];
        end
        
        %%% DEPENDS ON WHICH MODEL
        if strcmp(XRDIMAGE.Instr.dettype, '0')
            mapped_tth  = GeometricModelXRD0(...
                XRDIMAGE.Instr.centers./1000, ...
                XRDIMAGE.Instr.distance, ...
                XRDIMAGE.Instr.gammaY, XRDIMAGE.Instr.gammaX, ...
                Pixel2mm(polimg.radius', XRDIMAGE.Instr.pixelsize), polimg.azimuth, XRDIMAGE.Instr.detpars)';
        elseif strcmp(XRDIMAGE.Instr.dettype, '1')
            mapped_tth  = GeometricModelXRD1(...
                XRDIMAGE.Instr.centers./1000, ...
                XRDIMAGE.Instr.distance, ...
                XRDIMAGE.Instr.gammaY, XRDIMAGE.Instr.gammaX, ...
                Pixel2mm(polimg.radius', XRDIMAGE.Instr.pixelsize), polimg.azimuth, XRDIMAGE.Instr.detpars)';
        elseif strcmp(XRDIMAGE.Instr.dettype, '2')
            mapped_tth  = GeometricModelXRD2(...
                XRDIMAGE.Instr.centers./1000, ...
                XRDIMAGE.Instr.distance, ...
                XRDIMAGE.Instr.gammaY, XRDIMAGE.Instr.gammaX, ...
                Pixel2mm(polimg.radius', XRDIMAGE.Instr.pixelsize), polimg.azimuth, XRDIMAGE.Instr.detpars)';
        elseif strcmp(XRDIMAGE.Instr.dettype, '2a')
            mapped_tth  = GeometricModelXRD2a(...
                XRDIMAGE.Instr.centers./1000, ...
                XRDIMAGE.Instr.distance, ...
                XRDIMAGE.Instr.gammaY, XRDIMAGE.Instr.gammaX, ...
                Pixel2mm(polimg.radius', XRDIMAGE.Instr.pixelsize), polimg.azimuth, XRDIMAGE.Instr.detpars)';
        elseif strcmp(XRDIMAGE.Instr.dettype, '2b')
            mapped_tth  = GeometricModelXRD2b(...
                XRDIMAGE.Instr.centers./1000, ...
                XRDIMAGE.Instr.distance, ...
                XRDIMAGE.Instr.gammaY, XRDIMAGE.Instr.gammaX, ...
                Pixel2mm(polimg.radius', XRDIMAGE.Instr.pixelsize), polimg.azimuth, XRDIMAGE.Instr.detpars)';
        end
        
        polimg.mapped_tth_for_intensity = mapped_tth;
        
        [tth_grid, intensity_in_tth_grid]   = MapIntensityToTThGrid(XRDIMAGE, polimg);
        polimg.tth_grid                 = tth_grid;
        polimg.intensity_in_tth_grid    = intensity_in_tth_grid;
        
        figure(101)
        subplot(1,2,1)
        imagesc(log(abs(polimg.intensity_in_tth_grid))), axis square tight
        title('Caked image // radial position is corrected')
        hold off
        
        disp('###########################')
        disp(sprintf('mean pseudo-strain using p0 : %f', mean(abs(strain0(:)))))
        disp(sprintf('mean pseudo-strain using p  : %f\n', mean(abs(strain(:)))))
        disp(sprintf('std pseudo-strain using p0 : %f', std(strain0(:))))
        disp(sprintf('std pseudo-strain using p  : %f\n', std(strain(:))))
        
        %%% ASSIGN NEW INSTRUMENT PARAMETERS USING OPTIMIZATION RESULTS
        Instr   = XRDIMAGE.Instr;
        Instr.centers   = p(1:2);
        Instr.distance  = p(3);
        Instr.gammaX    = p(4);
        Instr.gammaY    = p(5);
        Instr.detpars   = p(6:end);
        
        if Analysis_Options.save_instrpars
            disp('###########################')
            disp(sprintf('Saving optimized innstrument parameters in %s\n', pfname_instr))
            save(pfname_instr, 'Instr')
        else
            disp('###########################')
            disp(sprintf('NOT saving optimized innstrument parameters for %s\n', pfname{i,1}))
        end
    end
end
