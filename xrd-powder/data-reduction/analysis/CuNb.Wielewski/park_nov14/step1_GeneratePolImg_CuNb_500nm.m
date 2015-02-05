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

% %%% MATERIAL PARAMETERS
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

%%% CALCULATE THEORETICAL TTH
% [d, th] = PlaneSpacings(XRDIMAGE.Material.lattparms, ...
%     'cubic', XRDIMAGE.Material.hkls', ...
%     XRDIMAGE.Instr.wavelength);
% tth     = 2*th;
% 
% XRDIMAGE.Material.tth       = tth;
% XRDIMAGE.Material.d_spacing = d;

%%% DATA REDUCTION FLAGS
Analysis_Options.make_polimg    = 1;
Analysis_Options.save_polimg    = 1;
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
        for j = 1:10:XRDIMAGE.Image.numframe
            disp('###########################')
            fprintf('Looking at %s frame number %d\n', pfname{i,1}, j)
            disp('###########################')
            
            pfname_polimage = [pfname{i,1}, '.frame', num2str(j), '.cor.polimg.mat'];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% POLAR REBINNING IF NECESSARY
            if XRDIMAGE.Image.corrected
                imgi    = ReadSUM(pfname{i,1});
            else
                imgi    = bg.*0;
                for k = 0:1:9
                    imgk    = NreadGE(pfname{i,1}, j + k);
                    imgi    = imgi + imgk;
                end
                imgi    = double(imgi - 10.*bg);
            end
            
            % NOTE ROTATION IS NEED TO ENSURE PROPER REBINNING
            imgi    = rot90(imgi);
            
            figure(1)
            imagesc(log(abs(imgi)))
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
end