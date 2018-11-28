function imdata = ReadDexela(filename, varargin)
% ReadDexela - reads .tiff images from Dexela.
%
%   USAGE:
%
%   imdata = ReadDexela(filename)
%
%   imdata = ReadDexela(filename, varargin)
%
%   INPUT:
%
%   filename            File name of the Dexela tiff file.
%
%   Input variable "filename" can be followed by a list of
%   parameter/value pairs which control certain plotting
%   features.  Options are:
%
%   Orientation         detector orientation. options are
%                       cable_UP - cable pointing along +Y (default)
%                       cable_DOWN - cable pointing along -Y
%                       cable_IB - cable pointing -X
%                       cable_OB - cable pointing +X
%
%   PlotImage           plot image data

% %%% SYNTHETIC IMAGE
% xxx = (1:1:2048)./204.8;
% yyy = (1:1:2048)./(204.8*2);
% imgi    = xxx'*yyy;
% imgi(1000:1200,:)   = 27;

% DEFAULT OPTIONS
optcell = {...
    'Orientation', 'cable_UP', ...
    'PlotImage', false, ...
    };

% UPDATE OPTION
opts    = OptArgs(optcell, varargin);

% READ IMAGE FILE
disp(sprintf('reading %s in %s orientation', filename, opts.Orientation));
imdata  = imread(filename);
if strcmpi(opts.Orientation, 'cable_UP')
    imdata  = rot90(imdata, -1);    %%% THIS MAKES WHAT YOU SEE WHAT YOU GET
elseif strcmpi(opts.Orientation, 'cable_DOWN')
    imdata  = rot90(imdata, 1);
elseif strcmpi(opts.Orientation, 'cable_IB')
    imdata  = rot90(imdata, 0);
elseif strcmpi(opts.Orientation, 'cable_OB')
    imdata  = rot90(imdata, 2);
else
    disp('no such detector orientation option')
    return
end

if opts.PlotImage
    imdata_disp = log(double(imdata));% cable_UP is as is
    
    xmin    = 0;
    ymin    = 0;
    [xmax, ymax] = size(imdata_disp);
    
    figure(1000)
    imagesc(rot90(imdata_disp,1)) %% DO NOT CHANGE
    colorbar horz
    caxis([6 7])
    colormap jet
    title({'loaded dexela image', 'check that the image matches the coordinate system'})
    xlabel('X_L (pixels)')
    ylabel('Y_L (pixels)')
    text(xmin, ymin, 'TO')
    text(xmax, ymin, 'TI')
    text(xmin, ymax, 'BO')
    text(xmax, ymax, 'BI')
    axis equal tight
    drawnow
end