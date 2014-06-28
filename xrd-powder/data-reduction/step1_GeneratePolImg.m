clear all
close all
clc

%%% INPUT PARAMETERS
XRDIMAGE.DarkField.pname    = '.\example\APS';
XRDIMAGE.DarkField.fbase    = 'LaB6_';
XRDIMAGE.DarkField.fnumber  = 4113;
XRDIMAGE.DarkField.numframe = 1;
XRDIMAGE.DarkField.numdigs  = 5;
XRDIMAGE.DarkField.fext     = 'ge2';

XRDIMAGE.Image.pname        = '.\example\APS';
XRDIMAGE.Image.fbase        = 'LaB6_';
XRDIMAGE.Image.fnumber      = [4116; 4117]; % 4116 / 4117
XRDIMAGE.Image.numframe     = 10;
XRDIMAGE.Image.numdigs      = 5;
XRDIMAGE.Image.fext         = 'ge2';

XRDIMAGE.Calib.pname        = '.\example\APS';
XRDIMAGE.Calib.fbase        = 'LaB6_';
XRDIMAGE.Calib.fnumber      = [4116; 4117]; % 4116 / 4117

%%% INSTRUMENT PARAMETERS GETS LOADED

%%% CAKE PARAMETERS
XRDIMAGE.CakePrms.bins(1)   = 72;           % number of azimuthal bins
XRDIMAGE.CakePrms.bins(2)   = 1500;         % number of radial bins
XRDIMAGE.CakePrms.bins(3)   = 20;           % number of angular bins
XRDIMAGE.CakePrms.origin(1) = 1023.628;         % x center in pixels, fit2d Y 
XRDIMAGE.CakePrms.origin(2) = 1019.544;         % y center in pixels, fit2d X
XRDIMAGE.CakePrms.sector(1) = -360/XRDIMAGE.CakePrms.bins(1)/2;     % start azimuth (min edge of bin) in degrees
XRDIMAGE.CakePrms.sector(2) = 360-360/XRDIMAGE.CakePrms.bins(1)/2;  % stop  azimuth (max edge of bin) in degrees
XRDIMAGE.CakePrms.sector(3) = 200;  % start radius (min edge of bin) in pixels
XRDIMAGE.CakePrms.sector(4) = 950;  % stop  radius (max edge of bin) in pixels
XRDIMAGE.CakePrms.azim      = 0:360/XRDIMAGE.CakePrms.bins(1):XRDIMAGE.CakePrms.sector(2);

%%% MATERIAL PARAMETERS
XRDIMAGE.Material.num       = 1;
XRDIMAGE.Material.lattparms = 4.1569162;        % LaB6
XRDIMAGE.Material.structure = 'simplecubic';
XRDIMAGE.Material.numpk     = 8;
XRDIMAGE.Material.pkrange    = [...
    2.7  3.8 4.7 6.1 6.7 8.2 8.7 9.1; ...
    2.95 4.1 5.0 6.7 7.0 8.5 9.0 9.4; ...
    ];
XRDIMAGE.Material.pkidx     = {...
    [1] [2] [3] [5] [6] [8] [9] [10] 
    };
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
pfname  = GenerateGEpfname(XRDIMAGE.DarkField);
bg      = NreadGE(pfname{1,1}, 1);

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
        
        imgi    = bg.*0;
        for j = 1:1:XRDIMAGE.Image.numframe
            imgj    = NreadGE(pfname{i,1}, j);
            imgi    = imgi + imgj;
        end
        imgi    = imgi - bg.*XRDIMAGE.Image.numframe;
        imgi    = rot90(imgi);
        
        % NOTE ROTATION IS NEED TO ENSURE PROPER REBINNING
        % img     = fliplr(rot90(img,1));
        
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