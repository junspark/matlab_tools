function polimg = PolarBinXRD(mesh, instr, cakeParms, img, varargin)
% PolarBinXRD - Polar integration of image data
%
%   USAGE:
%
%   polimg = PolarBinXRD(mesh, instr, cakeParms, img)
%
%   INPUT:
%
%   mesh
%       mesh structure for finite element integration
%
%   instr
%       instrument parameters to correct for the experimental geometry
%
%   cakeParms
%       caking parameters for integration
%
%   img
%       image data for integration
%
%   OUTPUT:
%
%   polimg
%       radially integrated data organized in struct variable
%
%   NOTE:
%
%   1) this does not do retangular pixels.

% default options
optcell = {...
    'PlotProgress', 'on', ...
    'DisplayProgress', 'on', ...
    };

% update option
opts    = OptArgs(optcell, varargin);

img     = double(img); 
imgi    = img;

Lx  = instr.numpixelsHorz;
Ly  = instr.numpixelsVert;

% !!! THESE ARE IN THE CARTESIAN FRAME !!!
x0  = cakeParms.origin(1);   % in pixels
y0  = cakeParms.origin(2);   % in pixels

x0plt   = x0;
y0plt   = instr.numpixelsVert - y0;

startAzi    = cakeParms.sector(1);
endAzi      = cakeParms.sector(2);

startRho    = cakeParms.sector(3);
endRho      = cakeParms.sector(4);

%%% NUMBER OF AZIMUTHAL BINS OVER ANGULAR RANGE DEFINED BY cakeParms.sector(1) AND cakeParms.sector(2)
numAzi  = cakeParms.bins(1);
%%% NUMBER OF RADIAL POINTS PER AZIMHUTHAL BIN OVER RADIAL RANGE DEFINED BY startRho AND endRho
numRho  = cakeParms.bins(2);
%%% NUMBER OF ETA POINTS PER AZIMUTH
numEta  = cakeParms.bins(3);

dAzi    = (endAzi - startAzi)/numAzi;
dEta    = dAzi/numEta;

R   = startRho:(endRho-startRho)/numRho:endRho;
R   = repmat(R, numEta + 1, 1)';

Rlist   = mean((R(1:end - 1,:) + R(2:end,:))./2,2);

polimg.azimuth   = cakeParms.azim;
polimg.radius    = zeros(numAzi, numRho);
polimg.intensity = zeros(numAzi, numRho);

if strcmpi(opts.PlotProgress, 'on')
    figure(1000)
    imagesc(log(abs(rot90(img,1))))
    hold on
    axis equal
    plot(x0plt, y0plt, 'rh')
    xlabel('X_L (pixels)')
    ylabel('Y_L (pixels)')
    drawnow
end
% keyboard
Xgrid   = 1:1:Lx;
Ygrid   = 1:1:Ly;
[Xgrid, Ygrid]  = meshgrid(Ygrid, Xgrid);

for ii = 1:1:numAzi
    if strcmpi(opts.DisplayProgress, 'on')
        fprintf('Processing sector %d of %d\n', ii, numAzi);
    end
    tic;
    
    azi_ini = polimg.azimuth(ii) - dAzi/2;
    azi_fin = polimg.azimuth(ii) + dAzi/2;
    
    TH  = azi_ini:dEta:azi_fin;
    TH  = repmat(TH, numRho + 1, 1);
    
    [x, y]	= pol2cart(deg2rad(TH),R);
    x = x0 + x; 
    y = y0 + y;
    
    if strcmpi(opts.PlotProgress, 'on')
        THplt   = azi_ini:dEta:azi_fin;
        THplt   = repmat(-THplt, numRho + 1, 1);
        
        [xplt, yplt]    = pol2cart(deg2rad(THplt),R);
        xplt    = x0plt + xplt;
        yplt    = y0plt + yplt;
        
        figure(1000)
        title(num2str(polimg.azimuth(ii)))
        plot(xplt, yplt, 'k.')
        drawnow
    end
    
    tic;
    if all(x(:) < Lx) && all((x(:) > 0)) && all(y(:) < Ly) && all(y(:) > 0)
        warn_user   = 0;
    else
        warn_user   = 1;
    end
    % keyboard
    V   = interp2(Xgrid, Ygrid, double(imgi), y, x, 'cubic');
    V   = mean(V,2);
    V   = (V(1:end-1) + V(2:end)) / 2;
    
    polimg.radius(ii,:)    = Rlist;
    polimg.intensity(ii,:) = V;
    
    dtime   = toc;
    if warn_user
        disp('Some requested nodal points are out of grid.')
    end
    if strcmpi(opts.DisplayProgress, 'on')
        fprintf('Processing time for sector %d is %1.4f\n', ii, dtime);
    end
end

if strcmpi(opts.PlotProgress, 'on')
    figure(1000)
    hold off
end
