clear all
close all
clc

%%% INPUT PARAMETERS
XRDIMAGE.DarkField.pname    = 'W:\__eval\CuNb_oct14';
XRDIMAGE.DarkField.fbase    = 'dark_0.3s_';
XRDIMAGE.DarkField.fnumber  = 28;
XRDIMAGE.DarkField.numframe = 1;
XRDIMAGE.DarkField.numdigs  = 5;
XRDIMAGE.DarkField.fext     = 'ge3';

XRDIMAGE.Image.pname        = 'W:\__eval\CuNb_oct14';
XRDIMAGE.Image.fbase        = 'CuNb_30nm_';
XRDIMAGE.Image.fnumber      = [351: 390]; % 50nm: 392 : 424 % 30 nm : 351: 390
XRDIMAGE.Image.numframe     = 180;
XRDIMAGE.Image.numdigs      = 5;
XRDIMAGE.Image.fext         = 'ge3';
XRDIMAGE.Image.corrected    = 0;

XRDIMAGE.Calib.pname        = 'W:\__eval\CuNb_oct14';
XRDIMAGE.Calib.fbase        = 'CeO2_1.5s_';
XRDIMAGE.Calib.fnumber      = [336; 337];
XRDIMAGE.Calib.fext         = 'ge3';

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
XRDIMAGE.CakePrms.bins(1)   = 24;           % number of azimuthal bins
XRDIMAGE.CakePrms.bins(2)   = 2000;         % number of radial bins
XRDIMAGE.CakePrms.bins(3)   = 1;            % number of angular bins
XRDIMAGE.CakePrms.origin(1) = 1024.190;         % x center in pixels, fit2d Y
XRDIMAGE.CakePrms.origin(2) = 1037.110;         % y center in pixels, fit2d X
XRDIMAGE.CakePrms.sector(1) = -360/XRDIMAGE.CakePrms.bins(1)/2;     % start azimuth (min edge of bin) in degrees
XRDIMAGE.CakePrms.sector(2) = 360-360/XRDIMAGE.CakePrms.bins(1)/2;  % stop  azimuth (max edge of bin) in degrees
XRDIMAGE.CakePrms.sector(3) = 200;  % start radius (min edge of bin) in pixels
XRDIMAGE.CakePrms.sector(4) = 950;  % stop  radius (max edge of bin) in pixels
XRDIMAGE.CakePrms.azim      = 0:360/XRDIMAGE.CakePrms.bins(1):XRDIMAGE.CakePrms.sector(2);
XRDIMAGE.CakePrms.fastint   = 1;

%%% MATERIAL PARAMETERS
XRDIMAGE.Material.pkbck     = 2;
XRDIMAGE.Material.pkfunc    = 4;
XRDIMAGE.Material.num           = 2;
XRDIMAGE.Material.lattparms{1}  = 3.6149;         % Cu
XRDIMAGE.Material.lattparms{2}  = 3.3004;         % Nb
XRDIMAGE.Material.structure{1}  = 'fcc';
XRDIMAGE.Material.structure{2}  = 'bcc';
XRDIMAGE.Material.hkls{1}       = load('fcc.hkls');
XRDIMAGE.Material.hkls{2}       = load('bcc.hkls');
XRDIMAGE.Material.latticetype{1}    = 'cubic';
XRDIMAGE.Material.latticetype{2}    = 'cubic';
for i = 1:1:XRDIMAGE.Material.num
    %%% CALCULATE THEORETICAL TTH
    [d, th] = PlaneSpacings(XRDIMAGE.Material.lattparms{i}, ...
        XRDIMAGE.Material.latticetype{i}, XRDIMAGE.Material.hkls{i}', ...
        XRDIMAGE.Instr.wavelength);
    tth	= 2*th;
    
    XRDIMAGE.Material.tth{i}        = tth;
    XRDIMAGE.Material.d_spacing{i}  = d;
end
XRDIMAGE.Material.pkrange       = [ ...
    4.206 4.996; ...
    5.059 5.359; ...
    5.777 6.276; ...
    6.391 6.791; ...
    7.797 8.321; ...
    8.313 8.713; ...
    9.039 9.618; ...
    9.825 10.15; ...
    ]';

for i = 1:1:XRDIMAGE.Material.num
    for j = 1:1:size(XRDIMAGE.Material.pkrange,2)
        idx1    = XRDIMAGE.Material.tth{i} > XRDIMAGE.Material.pkrange(1,j);
        idx2    = XRDIMAGE.Material.tth{i} < XRDIMAGE.Material.pkrange(2,j);
        idx     = find(idx1 & idx2);
        XRDIMAGE.Material.pkidx{i,j}    = idx;
    end
end

%%% DATA REDUCTION FLAGS
Analysis_Options.make_polimg    = 0;
Analysis_Options.save_polimg    = 0;
Analysis_Options.fits_spectra   = 1;
Analysis_Options.save_fits      = 1;
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
%%% LOAD XRD IMAGES
pfname  = GenerateGEpfname(XRDIMAGE.Image);
numimg  = length(pfname);
if Analysis_Options.fits_spectra
    for i = 1:1:numimg
        for j = 1:20:XRDIMAGE.Image.numframe
            disp('###########################')
            fprintf('Looking at %s frame number %d\n', pfname{i,1}, j)
            disp('###########################')
           
            pfname_polimage = [pfname{i,1}, '.frame', num2str(j), '.cor.polimg.mat'];
            pfname_pkfit    = [pfname{i,1}, '.frame', num2str(j), '.cor.pkfit.mat'];
            
            polimg  = load(pfname_polimage);
            polimg  = polimg.polimg;
            
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
            
            figure(2)
            subplot(1,2,1)
            imagesc(log(polimg.intensity_in_tth_grid )), axis square tight
            hold off
            
            subplot(1,2,2)
            plot(polimg.tth_grid, polimg.intensity_in_tth_grid)
            hold off
            
            disp('###########################')
            disp(sprintf('Fitting peaks in %s', pfname{i,1}))
            disp('###########################')
            
            x   = polimg.tth_grid(1,:);
            for k = 1:1:XRDIMAGE.CakePrms.bins(1)
                disp(sprintf('Looking at azimuthal bin %d of %d\n', k, XRDIMAGE.CakePrms.bins(1)))
                y   = polimg.intensity_in_tth_grid(k,:);
                
                figure(11)
                subplot(1,2,1)
                plot(x, y, 'k.')
                hold on
                for ii = 1:1:XRDIMAGE.Material.num
                    plot(XRDIMAGE.Material.tth{ii}, mean(y), '^')
                end
                axis([min(x) max(x) 0 max(y)+100])
                xlabel('2 theta (deg)')
                ylabel('intensity (arb. units)')
                title(['bin number : ', num2str(j)])
                
                for ii = 1:1:XRDIMAGE.Material.num
                    for jj = 1:1:size(XRDIMAGE.Material.pkrange,2)
                        if ~isempty(XRDIMAGE.Material.pkidx{ii,jj})
                            pkrange = XRDIMAGE.Material.pkrange(:,jj);
                            
                            idx = find(x >= pkrange(1) & x <= pkrange(2));
                            xr  = x(idx)';
                            yr  = y(idx)';
                            
                            pr0 = [...
                                max(yr)/10 ...
                                0.1 ...
                                0.5 ...
                                XRDIMAGE.Material.tth{ii}(XRDIMAGE.Material.pkidx{ii,jj}) ...
                                0 ...
                                3e2];
                            
                            lb  = [ ...
                                0 ...
                                0 ...
                                0 ...
                                XRDIMAGE.Material.tth{ii}(XRDIMAGE.Material.pkidx{ii,jj})-0.5 ...
                                -Inf ...
                                -Inf];
                            
                            ub  = [ ...
                                +Inf ...
                                +Inf ...
                                1 ...
                                XRDIMAGE.Material.tth{ii}(XRDIMAGE.Material.pkidx{ii,jj})+0.5 ...
                                +Inf ...
                                +Inf];
                            
                            y0  = pfunc(pr0,xr);
                            [pr, rsn, ~, ef]    = lsqcurvefit(@pfunc, pr0, xr, yr, ...
                                lb, ub, Analysis_Options.PkFitOptions);
                            yf  = pfunc(pr,xr);
                            
                            pkfit{ii}.hkl{XRDIMAGE.Material.pkidx{ii,jj},k}     = XRDIMAGE.Material.hkls{ii}(XRDIMAGE.Material.pkidx{ii,jj},:);
                            pkfit{ii}.amp(XRDIMAGE.Material.pkidx{ii,jj},k)     = pr(1);
                            pkfit{ii}.fwhm(XRDIMAGE.Material.pkidx{ii,jj},k)    = pr(2);
                            pkfit{ii}.mix(XRDIMAGE.Material.pkidx{ii,jj},k)     = pr(3);
                            pkfit{ii}.tth(XRDIMAGE.Material.pkidx{ii,jj},k)     = pr(4);
                            pkfit{ii}.bkg{XRDIMAGE.Material.pkidx{ii,jj},k}     = pr(5:end);
                            pkfit{ii}.rsn(XRDIMAGE.Material.pkidx{ii,jj},k)     = rsn;
                            pkfit{ii}.ef(XRDIMAGE.Material.pkidx{ii,jj},k)      = ef;
                            pkfit{ii}.rwp(XRDIMAGE.Material.pkidx{ii,jj},k)     = ErrorRwp(yr, yf);
                           
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
                            title(['peak number : ', num2str(XRDIMAGE.Material.pkidx{ii,jj}), ' of ', num2str(ii)])
                            hold off
                        end
                    end
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
end