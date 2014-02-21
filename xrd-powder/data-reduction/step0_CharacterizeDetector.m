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

%%% INSTRUMENT PARAMETERS
XRDIMAGE.Instr.energy       = 61.999;       % keV
XRDIMAGE.Instr.wavelength   = keV2Angstrom(XRDIMAGE.Instr.energy);  % wavelength (Angstrom)
XRDIMAGE.Instr.detectordim  = 409.6;        % mm
XRDIMAGE.Instr.pixelsize    = 0.2;          % mm
XRDIMAGE.Instr.distance     = 989.0390;     % mm
XRDIMAGE.Instr.centers      = [-986.7767 -21.6210]; % center offsets x & y
XRDIMAGE.Instr.gammaX       = -0.0036;
XRDIMAGE.Instr.gammaY       = -0.0047;
XRDIMAGE.Instr.dims         = XRDIMAGE.Instr.detectordim/XRDIMAGE.Instr.pixelsize;   % total number of rows in the full image

% RADIAL CORRECTION
% 0 : no correction
% 1 : constant radial offset
% 2 : PROPOSED BY ISSN 0909-0495 LEE
XRDIMAGE.Instr.dettype  = '0';

% 0 : []
% 1 : constant value
% 2 : [a1 a2 n1 n2 rhod]
XRDIMAGE.Instr.detpars  = [-3.174e-5 -2.595e-4 3.111 2.495 240.8 2];

%%% CAKE PARAMETERS
XRDIMAGE.CakePrms.bins(1)   = 72;           % number of azimuthal bins
XRDIMAGE.CakePrms.bins(2)   = 1500;         % number of radial bins
XRDIMAGE.CakePrms.bins(3)   = 45;           % number of angular bins
XRDIMAGE.CakePrms.origin(1) = 1024;         % x center in pixels, fit2d convention
XRDIMAGE.CakePrms.origin(2) = 1024;         % y center in pixels, fit2d convention
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

%%% DATA REDUCTION FLAGS
XRDIMAGE.Options.make_polimg    = 0;
XRDIMAGE.Options.save_polimg    = 0;
XRDIMAGE.Options.fits_spectra   = 0;
XRDIMAGE.Options.save_fits      = 0;
XRDIMAGE.Options.find_instrpars = 1;
XRDIMAGE.Options.save_instrpars = 0;
XRDIMAGE.Options.find_detpars	= 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOAD XRD IMAGES
%%% BACKGROUND
pfname      = GenerateGEpfname(XRDIMAGE.DarkField);
return
bg          = NreadGE(pfname, 1);

%%% XRD IMAGES
pfname

% pfname      = GenerateGEpfname(XRDIMAGE.Image);
numframe    = 1:1:XRDIMAGE.Image.numframe;
img         = bg.*0;
for i = numframe
    im  = NreadGE(pfname, i);
    % NEW GE USED AT APS 201103
    img = rot90(im-bg) + img;
end

figure(1)
imagesc(abs(log(img)))
colorbar
axis equal tight
xlabel('X_{L}')
ylabel('Y_{L}')
hold off
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% POLAR REBINNING
pf  = [pfname, '.polimg.mat'];
if XRDIMAGE.Options.make_polimg
    mesh  = BuildMeshXRD(XRDIMAGE.Instr.dims);
    
    % NOTE ROTATION IS NEED TO ENSURE PROPER REBINNING
    img     = fliplr(rot90(img,1));
    polimg  = PolarBinXRD(mesh, ...
        XRDIMAGE.Instr, ...
        XRDIMAGE.CakePrms, ...
        img);
    
    if XRDIMAGE.Options.save_polimg
        save(pf, 'polimg')
    end
else
    load(pf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PK FITTING
%%% FITTING OPTIONS
options = optimset(...
    'MaxIter', 5e5,...
    'MaxFunEvals',3e5);

%%% CALCULATE THEORETICAL TTH
pf      = [XRDIMAGE.Material.structure, '.hkls'];
hkls	= load(pf);
[~, th] = PlaneSpacings(XRDIMAGE.Material.lattparms, ...
    'cubic', hkls', ...
    XRDIMAGE.Instr.wavelength);
tth     = 2*th;

pf  = [pfname, '.polimg.fit.mat'];
if XRDIMAGE.Options.fits_spectra
    for i = 1:1:XRDIMAGE.CakePrms.bins(1)
        x   = polimg.radius(i,:);
        x   = Pixel2mm(x, XRDIMAGE.Instr.pixelsize);  % CONVERT TO MM FROM PIXELS
        y   = polimg.intensity(i,:);
        
        figure(10)
        plot(x,y,'k.-')
        hold on
        plot(XRDIMAGE.Instr.distance*tand(tth), 1000, 'g^')
        axis([min(x) max(x) 0 max(y)+100])
        
        for j = 1:1:XRDIMAGE.Material.numpk
            if i == 1
                pkrange = XRDIMAGE.Material.pkrange(:,j);
                pkrange = XRDIMAGE.Instr.distance.*tand(pkrange);
                idx = find(x >= pkrange(1) & x <= pkrange(2));
                xr  = x(idx)';
                yr  = y(idx)';
                
                pr0 = [...
                    max(yr)/5 ...
                    0.5 ...
                    0.15 ...
                    XRDIMAGE.Instr.distance*tand(tth(XRDIMAGE.Material.pkidx{j})) + 1 ...
                    0 ...
                    2e3];
            else
                pkrange = [pkfit.rho(i-1,j)-2.5 pkfit.rho(i-1,j)+2.5];
                idx = find(x >= pkrange(1) & x <= pkrange(2));
                xr  = x(idx)';
                yr  = y(idx)';
                
                pr0 = [...
                    pkfit.amp(i-1,j) ...
                    pkfit.mix(i-1,j) ...
                    pkfit.fwhm(i-1,j) ...
                    pkfit.rho(i-1,j) ...
                    pkfit.bkg{i-1,j}];
            end
            plot(xr, yr, 'b.')
            
            y0  = pfunc(pr0,xr);
            plot(xr, y0, 'r-')
            
            [pr, rsn, ~, ef]    = lsqcurvefit(@pfunc, pr0, xr, yr, ...
                [], [], options);
            yf  = pfunc(pr,xr);
            plot(xr, yf, 'g-')
            
            pkfit.amp(i,j)  = pr(1);
            pkfit.mix(i,j)  = pr(2);
            pkfit.fwhm(i,j) = pr(3);
            pkfit.rho(i,j)  = pr(4);
            pkfit.bkg{i,j}  = pr(5:end);
            pkfit.rsn(i,j)  = rsn;
            pkfit.ef(i,j)   = ef;
            pkfit.rwp(i,j)  = ErrorRwp(yr, yf);
            disp([i j])
        end
        hold off
    end
    if XRDIMAGE.Options.save_fits
        save(pf, 'pkfit')
    end
else
    load(pf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% APPLY/FIND GEOMETRICAL MODEL
%%% FITTING OPTIONS
if XRDIMAGE.Options.find_instrpars
    % 'TolFun', 1e-18, ...
    % 'TolX', 1e-22, ...
    options = optimset(...
        'DerivativeCheck', 'off', ...
        'MaxIter', 1e5, ...
        'MaxFunEvals',3e5, ...
        'TypicalX',[10 -10 1 1 1 XRDIMAGE.Instr.detpars], ...
        'Display','final');
    
    p0  = [...
        XRDIMAGE.Instr.centers, ...
        XRDIMAGE.Instr.distance, ...
        XRDIMAGE.Instr.gammaX, ...
        XRDIMAGE.Instr.gammaY, ...
        XRDIMAGE.Instr.detpars];
    p           = lsqnonlin(@ApplyGeometricModel,p0,[],[],options);
    
    if XRDIMAGE.Options.save_instrpars
        pf  = [pfname, '.instr.mat'];
        
        XRDIMAGE.Instr.centers  = p(1:2);
        XRDIMAGE.Instr.distance = p(3);
        XRDIMAGE.Instr.gammaX   = p(4);
        XRDIMAGE.Instr.gammaY   = p(5);
        XRDIMAGE.Instr.detpars  = p(6:end);
        
        Instr   = XRDIMAGE.Instr;
        save(pf, 'Instr');
    end
else
    pf  = GenerateGEpfname(XRDIMAGE.Calib,1);
    pf  = [pf, '.instr.mat'];
    
    load(pf)
    XRDIMAGE.Instr  = Instr;
    
    p   = [...
        XRDIMAGE.Instr.centers, ...
        XRDIMAGE.Instr.distance, ...
        XRDIMAGE.Instr.gammaX, ...
        XRDIMAGE.Instr.gammaY,  ...
        XRDIMAGE.Instr.detpars];
end

%%% PLOT INSTRUMENT PRMS RESULT
figure(100)
tth_meas    = ApplyGeometricModel(p) + ...
        repmat(tth([XRDIMAGE.Material.pkidx{:}]),XRDIMAGE.CakePrms.bins(1),1)';
    
tth_calc    = tth([XRDIMAGE.Material.pkidx{:}])';
tth_calc    = repmat(tth_calc, 1, XRDIMAGE.CakePrms.bins(1));
dth         = (tth_meas - tth_calc)./tth_calc;
for i = 1:1:XRDIMAGE.CakePrms.bins(1)
    plot(dth(:,i), 'k.')
    hold on
end
xlabel('peak number')
ylabel('d\theta/\theta_{calc}')
grid on
axis tight

figure(101)
lgnd_shape  = {'s', 'o', '^', 'v', 's', 'o', '^', 'v'};
lgnd_edge   = {'r', 'g', 'b', 'k', 'r', 'g', 'b', 'k'};
lgnd_face   = {'r', 'g', 'b', 'k', 'none', 'none', 'none', 'none'};
lgnd        = num2str([1:8]');
for i = 1:1:XRDIMAGE.Material.numpk
    plot(dth(i,:) + i*1e-4, XRDIMAGE.CakePrms.azim, ...
        'Marker', lgnd_shape{i}, ...
        'MarkerEdgeColor', lgnd_edge{i}, ...
        'MarkerFaceColor', lgnd_face{i}, ...
        'Color', lgnd_edge{i})
    hold on
end
legend(lgnd)
xlabel('d\theta/\theta_{calc}')
ylabel('\eta (deg)')
grid on
axis tight