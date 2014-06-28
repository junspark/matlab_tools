clear all
close all
clc

%%% INPUT PARAMETERS
XRDIMAGE.Image.pname        = '.\ge3\corrected';
XRDIMAGE.Image.fbase        = 'DiaPac_';
XRDIMAGE.Image.fnumber      = 80009:80050; % 4116 / 4117
XRDIMAGE.Image.numframe     = 1;
XRDIMAGE.Image.numdigs      = 5;
XRDIMAGE.Image.fext         = 'ge3.sum';

XRDIMAGE.Calib.pname        = '.\ge3\corrected';
XRDIMAGE.Calib.fbase        = 'DiaPac_';
XRDIMAGE.Calib.fnumber      = 80113; % 4116 / 4117
XRDIMAGE.Calib.fext         = 'ge3.sum';

%%% INSTRUMENT PARAMETERS GETS LOADED
Instr.energy        = 0;
Instr.wavelength	= 0;
Instr.detectorsize	= 0;
Instr.pixelsize     = 0;
Instr.distance      = 0;
Instr.centers       = [0, 0];
Instr.gammaX        = 0;
Instr.gammaY        = 0;
Instr.numpixels     = 0;
Instr.dettype       = '2a';
Instr.detpars       = [0 0 0 0 0 0];
pfname  = GenerateGEpfname(XRDIMAGE.Calib);
for i = 1:1:length(XRDIMAGE.Calib.fnumber)
    pfname_instr    = [pfname{i,1}, '.instr.mat'];
    Instri  = load(pfname_instr);
    Instri  = Instri.Instr;
    
    Instr.energy        = Instr.energy + Instri.energy;
    Instr.wavelength    = Instr.wavelength + Instri.wavelength;
    Instr.detectorsize  = Instr.detectorsize + Instri.detectorsize;
    Instr.pixelsize     = Instr.pixelsize + Instri.pixelsize;
    Instr.distance      = Instr.distance + Instri.distance;
    Instr.centers       = Instr.centers + Instri.centers;
    Instr.gammaX        = Instr.gammaX + Instri.gammaX;
    Instr.gammaY        = Instr.gammaY + Instri.gammaY;
    Instr.numpixels     = Instr.numpixels + Instri.numpixels;
    Instr.detpars       = Instr.detpars + Instri.detpars;
end
Instr.energy        = Instr.energy./length(XRDIMAGE.Calib.fnumber);
Instr.wavelength    = Instr.wavelength./length(XRDIMAGE.Calib.fnumber);
Instr.detectorsize  = Instr.detectorsize./length(XRDIMAGE.Calib.fnumber);
Instr.pixelsize     = Instr.pixelsize./length(XRDIMAGE.Calib.fnumber);
Instr.distance      = Instr.distance./length(XRDIMAGE.Calib.fnumber);
Instr.centers       = Instr.centers./length(XRDIMAGE.Calib.fnumber);
Instr.gammaX        = Instr.gammaX./length(XRDIMAGE.Calib.fnumber);
Instr.gammaY        = Instr.gammaY./length(XRDIMAGE.Calib.fnumber);
Instr.numpixels     = Instr.numpixels./length(XRDIMAGE.Calib.fnumber);
Instr.detpars       = Instr.detpars./length(XRDIMAGE.Calib.fnumber);

Instr.omega = 0;
Instr.chi   = 0;

XRDIMAGE.Instr  = Instr;

%%% CAKE PARAMETERS
XRDIMAGE.CakePrms.bins(1)   = 72;           % number of azimuthal bins
XRDIMAGE.CakePrms.bins(2)   = 1500;         % number of radial bins
XRDIMAGE.CakePrms.bins(3)   = 20;           % number of angular bins
XRDIMAGE.CakePrms.origin(1) = 1021.744;         % x center in pixels, fit2d Y 
XRDIMAGE.CakePrms.origin(2) = 987.765;         % y center in pixels, fit2d X
XRDIMAGE.CakePrms.sector(1) = -360/XRDIMAGE.CakePrms.bins(1)/2;     % start azimuth (min edge of bin) in degrees
XRDIMAGE.CakePrms.sector(2) = 360-360/XRDIMAGE.CakePrms.bins(1)/2;  % stop  azimuth (max edge of bin) in degrees
XRDIMAGE.CakePrms.sector(3) = 200;  % start radius (min edge of bin) in pixels
XRDIMAGE.CakePrms.sector(4) = 950;  % stop  radius (max edge of bin) in pixels
XRDIMAGE.CakePrms.azim      = 0:360/XRDIMAGE.CakePrms.bins(1):XRDIMAGE.CakePrms.sector(2);

%%% MATERIAL PARAMETERS
XRDIMAGE.Material.num       = 1;
XRDIMAGE.Material.lattparms = 3.5668;        % CeO2
XRDIMAGE.Material.structure = 'fcc';
XRDIMAGE.Material.numpk     = 5;
XRDIMAGE.Material.pkidx     = {...
    [1] [3] [4] [5] [6]
    };
XRDIMAGE.Material.pkrange    = [...
    4.8295    5.5927    7.9541    9.3471    9.7682   11.2996;
    5.0295    5.7927    8.1541    9.5471    9.9682   11.4996;
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

%%% DATA REDUCTION FLAGS
Analysis_Options.make_polimg    = 1;
Analysis_Options.save_polimg    = 1;
Analysis_Options.generateESG	= 1;
Analysis_Options.fits_spectra   = 0;
Analysis_Options.save_fits      = 0;
Analysis_Options.find_instrpars = 0;
Analysis_Options.save_instrpars = 0;
Analysis_Options.find_detpars	= 0;

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
pfname  = GenerateGEpfname(XRDIMAGE.Image);
numimg  = length(pfname);
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