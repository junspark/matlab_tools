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
%   SPECSettingCorrect  detector initialization in SPEC is correct
%                       (default = true). If true, Orientation parameter is
%                       ignored. If false, Orientation parameter is used to
%
%   Orientation         detector orientation. options are
%                       cable_UP - cable pointing along +Y
%                       cable_DOWN - cable pointing along -Y (default)
%                       cable_IB - cable pointing -X
%                       cable_OB - cable pointing +X
%
%   PlotImage           plot image data

% DEFAULT OPTIONS
optcell = {...
    'SPECSettingCorrect', true, ...
    'Orientation', 'cable_DOWN', ...
    'PlotImage', false, ...
    };

% UPDATE OPTION
opts    = OptArgs(optcell, varargin);

% READ IMAGE FILE
disp(sprintf('reading %s in %s orientation', filename, opts.Orientation));
imdata  = imread(filename);

if opts.SPECSettingCorrect
    disp(sprintf('specified det orientation %s ignored', opts.Orientation));
    imdata  = rot90(imdata, -1);
else
    %%% IF SPECSettingCorrect FALSE, THIS MAKES WHAT YOU SEE WHAT YOU GET
    %%% THIS NEEDS TO BE CHECKED
	disp(sprintf('Dexela is in %s configuration', opts.Orientation))
    if strcmpi(opts.Orientation, 'cable_DOWN')
        imdata  = rot90(imdata, -1);    
    elseif strcmpi(opts.Orientation, 'cable_UP')
        disp('WARNING!!!! CABLE UP = 1')
        imdata  = rot90(imdata, 1);
    elseif strcmpi(opts.Orientation, 'cable_IB')
        disp('WARNING!!!! still needs to be implemented')
        imdata  = rot90(imdata, 0);
    elseif strcmpi(opts.Orientation, 'cable_OB')
        disp('WARNING!!!! still needs to be implemented')
        imdata  = rot90(imdata, 2);
    else
        disp('no such detector orientation option')
        return
    end
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
