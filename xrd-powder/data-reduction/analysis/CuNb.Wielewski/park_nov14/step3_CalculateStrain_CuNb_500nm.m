clear all
close all
clc

%%% INPUT PARAMETERS
XRDIMAGE.DarkField.pname    = 'W:\__eval\park_nov14\ge';
XRDIMAGE.DarkField.fbase    = 'CuNb_500nm_';
XRDIMAGE.DarkField.fnumber  = 9026;
XRDIMAGE.DarkField.numframe = 1;
XRDIMAGE.DarkField.numdigs  = 5;
XRDIMAGE.DarkField.fext     = 'ge3';

XRDIMAGE.Image.pname        = 'W:\__eval\park_nov14\ge';
XRDIMAGE.Image.fbase        = 'CuNb_500nm_';
XRDIMAGE.Image.fnumber      = [9019:9022, 9024:9025]; % 50nm: 392 : 424 % 30 nm : 351; 390
XRDIMAGE.Image.numframe     = 180;
XRDIMAGE.Image.numdigs      = 5;
XRDIMAGE.Image.fext         = 'ge3';
XRDIMAGE.Image.corrected    = 0;

XRDIMAGE.Calib.pname        = 'O:\park_nov14';
XRDIMAGE.Calib.fbase        = 'CeO2_';
XRDIMAGE.Calib.fnumber      = [746; 747];
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
XRDIMAGE.CakePrms.bins(1)   = 36;           % number of azimuthal bins
XRDIMAGE.CakePrms.bins(2)   = 2000;         % number of radial bins
XRDIMAGE.CakePrms.bins(3)   = 20;           % number of angular bins
XRDIMAGE.CakePrms.origin(1) = 1038.159;         % x center in pixels, fit2d Y 
XRDIMAGE.CakePrms.origin(2) = 1013.426;         % y center in pixels, fit2d X
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
    3.228 3.653; ...
    3.748 4.020; ...
    4.336 4.595; ...
    4.801 4.986; ...
    5.847 6.119; ...
    6.391 6.791; ...
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
            title(num2str(i))
            hold off
            
            subplot(1,2,2)
            plot(polimg.tth_grid, polimg.intensity_in_tth_grid)
            hold off
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            pfname_pkfit    = [pfname{i,1}, '.frame', num2str(j), '.cor.pkfit.mat'];
            
            pkfit   = load(pfname_pkfit);
            pkfit   = pkfit.pkfit;
            
            for k = 1:1:XRDIMAGE.CakePrms.bins(1)
                disp(sprintf('Looking at azimuthal bin %d of %d\n', k, XRDIMAGE.CakePrms.bins(1)))
                
                for ii = 1:1:XRDIMAGE.Material.num
                    for jj = 1:1:size(XRDIMAGE.Material.pkrange,2)
                        if ~isempty(XRDIMAGE.Material.pkidx{ii,jj})
%                             pkfit{ii}.hkl{XRDIMAGE.Material.pkidx{ii,jj},k}
%                             pkfit{ii}.amp(XRDIMAGE.Material.pkidx{ii,jj},k)
%                             pkfit{ii}.tth(XRDIMAGE.Material.pkidx{ii,jj},k)
                            
                            th0 = XRDIMAGE.Material.tth{ii}(XRDIMAGE.Material.pkidx{ii,jj})/2;
                            th  = pkfit{ii}.tth(XRDIMAGE.Material.pkidx{ii,jj},k)/2;
                            eqq = sind(th) / sind(th0) - 1;
                            LatticeStrain{ii}(XRDIMAGE.Material.pkidx{ii,jj},k) = eqq;
                        end
                    end
                end
            end
            
            for ii = 1:1:XRDIMAGE.Material.num
                for jj = 1:1:size(XRDIMAGE.Material.pkrange,2)
                    if ~isempty(XRDIMAGE.Material.pkidx{ii,jj})
                        Amp = pkfit{ii}.amp(XRDIMAGE.Material.pkidx{ii,jj},:);
                        mu  = mean(Amp); sigma  = std(Amp);
                        AmpThresh   = mu - 0.0*sigma;
                        idxAmp  = Amp < 20;
                        
                        ef  = pkfit{ii}.ef(XRDIMAGE.Material.pkidx{ii,jj},:);
                        
                        Rwp = pkfit{ii}.rwp(XRDIMAGE.Material.pkidx{ii,jj},:);
                        idxRwp  = Rwp < 0.031;
                        
                        %%% NOTE THERE WAS A BUG IN INDEX SO SAVED FILES
                        %%% HAVE SWAPPED MIX AND FWHM PARAMETERS
                        mix     = pkfit{ii}.fwhm(XRDIMAGE.Material.pkidx{ii,jj},:);
                        fwhm    = pkfit{ii}.mix(XRDIMAGE.Material.pkidx{ii,jj},:);
                        idxfwhm = fwhm < 0.1;
                        
                        idx = idxAmp & idxRwp & idxfwhm;
                        idx = idxRwp;
                        
                        if ii == 1
                            MatStr  = 'Cu';
                        elseif ii == 2
                            MatStr  = 'Nb';
                        end
                        
                        figure(10*ii + jj)
                        subplot(2,2,1)
                        plot(XRDIMAGE.CakePrms.azim(:), Amp, 'o')
                        title(sprintf('%s \\{%d %d %d\\}', MatStr, pkfit{ii}.hkl{XRDIMAGE.Material.pkidx{ii,jj},1}))
                        axis tight
                        
                        subplot(2,2,2)
                        plot(XRDIMAGE.CakePrms.azim(:), Rwp, 'o')
                        title(sprintf('%s \\{%d %d %d\\}', MatStr, pkfit{ii}.hkl{XRDIMAGE.Material.pkidx{ii,jj},1}))
                        axis tight
                        
                        subplot(2,2,3)
                        scatter(XRDIMAGE.CakePrms.azim(:), LatticeStrain{ii}(XRDIMAGE.Material.pkidx{ii,jj},:), ...
                            30, Amp(:))
                        title(sprintf('%s \\{%d %d %d\\}', MatStr, pkfit{ii}.hkl{XRDIMAGE.Material.pkidx{ii,jj},1}))
                        hold on
                        axis([0 360 -2.5e-3 2.5e-3])
                        
                        subplot(2,2,4)
                        scatter(XRDIMAGE.CakePrms.azim(:), LatticeStrain{ii}(XRDIMAGE.Material.pkidx{ii,jj},:), ...
                            30, Rwp(:))
                        title(sprintf('%s \\{%d %d %d\\}', MatStr, pkfit{ii}.hkl{XRDIMAGE.Material.pkidx{ii,jj},1}))
                        hold on
                        axis([0 360 -2.5e-3 2.5e-3])
                        
%                         figure(100*ii + jj)
%                         subplot(2,2,1)
%                         plot(XRDIMAGE.CakePrms.azim(idx), Amp(idx), 'o')
%                         title(sprintf('%s \\{%d %d %d\\}', MatStr, pkfit{ii}.hkl{XRDIMAGE.Material.pkidx{ii,jj},1}))
%                         axis tight
%                         
%                         subplot(2,2,2)
%                         plot(XRDIMAGE.CakePrms.azim(idx), Rwp(idx), 'o')
%                         title(sprintf('%s \\{%d %d %d\\}', MatStr, pkfit{ii}.hkl{XRDIMAGE.Material.pkidx{ii,jj},1}))
%                         axis tight
%                         
%                         subplot(2,2,3)
%                         scatter(XRDIMAGE.CakePrms.azim(idx), LatticeStrain{ii}(XRDIMAGE.Material.pkidx{ii,jj},idx), ...
%                             30, Amp(idx))
%                         title(sprintf('%s \\{%d %d %d\\}', MatStr, pkfit{ii}.hkl{XRDIMAGE.Material.pkidx{ii,jj},1}))
%                         hold on
%                         axis([0 360 -1e-3 1e-3])
%                         
%                         subplot(2,2,4)
%                         scatter(XRDIMAGE.CakePrms.azim(idx), LatticeStrain{ii}(XRDIMAGE.Material.pkidx{ii,jj},idx), ...
%                             30, Rwp(idx))
%                         title(sprintf('%s \\{%d %d %d\\}', MatStr, pkfit{ii}.hkl{XRDIMAGE.Material.pkidx{ii,jj},1}))
%                         hold on
%                         axis([0 360 -1e-3 1e-3])
                    end
                end
            end
%             pause
% return
        end
    end
end