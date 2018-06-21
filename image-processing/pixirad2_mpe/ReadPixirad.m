function [csq, pixelsize] = ReadPixirad(pfname, varargin)
% ReadPixirad - read Pixirad hexagonal grid file.
%
%   INPUT:
%
%   pfname
%       name of the Pixirad image file.
%
%   version
%       pixirad version. 'pixi1' is for 1 panel pixirad from the APS
%       detector pool. 'pixi2' is for 2 panel pixirad at APS 1-ID beamline.
%       In the case of 'pixi1', it is assumed that the images are saved with
%       correction already applied. In the case of 'pixi2', it is assumed
%       that the images are saved without the correction and the correction
%       is applied in this code. Only the 1 color mode correction is
%       implemented. It is also assumed that the images are acquired with
%       AD Transform plug in with Rot270 option enabled. 'pixi2' is for any
%       images taken with the 1-id pixirad2 detector before mark rivers
%       fixed the driver in 2018-1. all images takens after 2018-1 should
%       use 'pixi3' option (which is same as pixi2 without the arbirary
%       gaussian filtering).
%
%   nc (optional - pixi1)
%       number of horizontal nodes (default = 476).
%
%   nr (optional - pixi1)
%       number of vertical nodes (default = 512).
%
%   nxsq (optional - pixi1)
%       number of pixels along x in the output square grid data. number of
%       pixels along y is computed based on this number to make the pixel
%       square.
%
%   display (optional - pixi1)
%       displays the sqaure grid image for confirmation (default = off).
%
%   apply_ct (optional - pixi2)
%       applies correction table if 1; default value is 0 (no).
%
%   pfname_ct (optional - pixi2)
%       full path and file name of the correction table. If not provided,
%       correction table located in 
%       /home/beams/S1IDUSER/mnt/s1b/__eval/matlab_tools/image-processing/pixirad2_mpe/sensitivity_map_Pb.mat
%       This file was generated from pixirad_jul17 data using method
%       described in pixirad_jul17_summary.pdf
%
%   apply_bpt (optional - pixi2)
%       applies bad pixel table if 1; default value is 0 (no).
%
%   pfname_bpt (optional - pixi2)
%       full path and file name of the bad pixel table. If not provided,
%       correction table located in 
%       /home/beams/S1IDUSER/mnt/s1b/__eval/matlab_tools/image-processing/pixirad2_mpe/bad_pixel_map_Pb.mat
%       This file was generated from pixirad_jul17 data using method
%       described in pixirad_jul17_summary.pdf
%       Pixels labeled as bad are visited iteratively (as there are many
%       bad pixels next to each other ... streaks) until they are all
%       visited at least once.
%
%   apply_gaussian_filter (optional - pixi2)
%       apply gaussian filter to remove spikes and noise; default value is 0 (no), ...
%
%   gaussian_filter_sigma (optional - pixi2)
%       sigma of the gaussian filter to be applied; default value is 2 and
%       should not be too large
%
%   OUTPUT:
%
%   csq
%       Pixirad image file information mapped to an image with square
%       pixels
%
%   pixelsize
%       Pixel size of the square pixels in mm. When nxsq is 476, resulting
%       pixel size of the square pixel is approximately 0.052 mm.
%       Ultimately, this needs to be found from optimization.

% default options
optcell = {...
    'version', 'pixi1', ...
    'nc', 476, ...
    'nr', 512, ...
    'nxsq', 476, ...
    'display', 'off', ...
    'apply_ct', 0, ...
    'pfname_ct', '/home/beams/S1IDUSER/mnt/s1b/__eval/matlab_tools/image-processing/pixirad2_mpe/sensitivity_map_Pb.mat', ...
    'apply_bpt', 0, ...
    'pfname_bpt', '/home/beams/S1IDUSER/mnt/s1b/__eval/matlab_tools/image-processing/pixirad2_mpe/bad_pixel_map_Pb.mat', ...
    'apply_gaussian_filter', 0, ...
    'gaussian_filter_sigma', 2, ...
    };

% update option
opts    = OptArgs(optcell, varargin);

if strcmpi(opts.version, 'pixi1')
    % read in image
    imdata      = double(imread(pfname));
    [nr, nc]    = size(imdata);
    
    if nc ~= opts.nc
        disp(sprintf('user input or default : %d', opts.nc))
        disp(sprintf('image size in x       : %d', nc))
        error('number of columns does not match')
    elseif nr ~= opts.nr
        disp(sprintf('user input or default : %d', opts.nr))
        disp(sprintf('image size in y       : %d', nr))
        error('number of rows does not match')
    else
        disp(sprintf('%d x %d image', nc, nr))
    end
    
    % WHEN nxsq == 476; pixelsize is 0.052 mm per pixel
    pixelsize   = opts.nxsq/476*0.052;
    
    pfname_xymap    = ['pixirad.map.nc.', num2str(nc), '.nr.', num2str(nr), '.mat'];
    if exist(pfname_xymap, 'file')
        load(pfname_xymap)
    else
        disp(sprintf('pixirad map %s does not exist!', pfname_xymap))
        disp(sprintf('creating %s.', pfname_xymap))
        x1  = 1:1:nc;
        x2  = 0.5:1:(nc-0.5);
        
        x   = [];
        y   = [];
        
        ct  = 1;
        for i = 1:1:nr
            if mod(i,2) == 1
                x   = [x; x1'];
                y   = [y; ct*ones(length(x1),1)./sind(60)];
            else
                x   = [x; x2'];
                y   = [y; ct*ones(length(x2),1)./sind(60)];
            end
            ct  = ct + 1;
        end
        save(pfname_xymap, 'x', 'y')
    end
    
    imdata  = imdata';
    imdata  = imdata(:);
    
    %%% CREATE INTERPOLATION
    F   = TriScatteredInterp(x, y, imdata);
    
    %%% GENERATE SQUARE GRID
    dx  = (max(x) - min(x))/(opts.nxsq - 1);
    dy  = dx;
    xsq = min(x):dx:max(x);
    ysq = min(y):dy:max(y);
    
    [xsq, ysq]  = meshgrid(xsq, ysq);
    
    %%% MAP HEX GRID DATA TO SQUARE GRID
    csq = F(xsq, ysq);
    
    if strcmpi(opts.display, 'on')
        figure(1000)
        subplot(1,3,1)
        imagesc(log(imdata))
        axis equal
        
        subplot(1,3,2)
        scatter(x, y, 10, log(imdata))
        axis equal tight
        
        subplot(1,3,3)
        scatter(xsq(:), ysq(:), 10, log(csq(:)))
        axis equal tight
    end
elseif strcmpi(opts.version, 'pixi2')
    % read in image
    csq = double(imread(pfname));
    
    % ONLY CRRM MODE (PIXEL MODE) SUPPORTED - THIS IS GENERATED BASED ON pixirad_jul17 DATA
    % ct  = load(opts.pfname_ct);
    % ct  = ct.correction_map;
    % ct  = fliplr(ct);
    
    % ONLY CRRM MODE (PIXEL MODE) SUPPORTED - THIS IS FROM THE VENDOR
    % ct  = imread('/home/beams/S1IDUSER/mnt/s1a/misc/pixirad2/usb.after_repair/Calibrations/2010_crrm.tif');
    % ct  = fliplr(ct');
    
    ct  = load(opts.pfname_ct);
    ct  = ct.ct;
    
    bpt = load(opts.pfname_bpt);
    bpt = bpt.bpt;
    
    %%% SENSITIVITY CORRECTOIN
    if opts.apply_ct
        % csq = csq./ct;
        csq = csq.*ct;
    end
    
    %%% BAD PIXEL CORRECTION
    %%% MULTIPLE PASSES TO REMOVE ALL THE BAD PIXELS
    %%% MIGHT BE PATH DEPENDENT
    if opts.apply_bpt
        % [n, m]  = size(csq);
        
        there_are_bps   = true;
        counter = 0;
        while there_are_bps && (counter < 5)
            disp(sprintf('removing bad pixels - pass number %d', counter))
            disp(sprintf('removing bad pixels - %d bad pixels exist', sum(bpt(:))))
            [x_bpt, y_bpt]      = find(bpt == 1);
            
            %%% GET NEIGHBOR PIXEL VALUES
            for i = 1:1:length(x_bpt)
%                 pv  = nan(8,1);
%                 if (x_bpt(i)-1) >= 1
%                     if (y_bpt(i)-1) >= 1
%                         pv(1)   = csq(x_bpt(i)-1, y_bpt(i)-1); % nw
%                     end
%                     if (y_bpt(i)+1) <= m
%                         pv(2)   = csq(x_bpt(i)-1, y_bpt(i)+1); % ne
%                     end
%                     pv(3)   = csq(x_bpt(i)-1, y_bpt(i)); % n
%                 end
%                 
%                 if (x_bpt(i)+1) <= n
%                     if (y_bpt(i)-1) >= 1
%                         pv(4)   = csq(x_bpt(i)+1, y_bpt(i)-1); % sw
%                     end
%                     if (y_bpt(i)+1) <= m
%                         pv(5)   = csq(x_bpt(i)+1, y_bpt(i)+1); % se
%                     end
%                     pv(6)   = csq(x_bpt(i)+1, y_bpt(i)); % s
%                 end
%                 
%                 if (y_bpt(i)+1) <= m
%                     pv(7)   = csq(x_bpt(i), y_bpt(i)+1); % e
%                 end
%                 if (y_bpt(i)-1) >= 1
%                     pv(8)   = csq(x_bpt(i), y_bpt(i)-1); % w
%                 end
%                 pv(find(pv == 0))   = nan;
                % csq(x_bpt(i), y_bpt(i)) = mean(pv, 'omitnan');
                
                imc = imcrop(csq, [y_bpt(i)-1 x_bpt(i)-1 2 2]);
                imc = imc(:);
                imc(find(imc == 0)) = nan;
                
                %%% SOME USERS DO NOT HAVE OMITNAN
                if verLessThan('matlab', '9.0')
                    disp('this will work for now but upgrade matlab')
                    idx = ~isnan(imc);
                    csq(x_bpt(i), y_bpt(i)) = sum(imc(idx))/sum(idx);
                else
                    csq(x_bpt(i), y_bpt(i)) = mean(imc, 'omitnan');
                end
            end
            bpt = (csq == 0);
            
            there_are_bps = ~isempty(find(bpt == 1));
            
            counter = counter + 1;
        end
    end
    
    %%% APPLY GAUSSIAN FILTER
    if opts.apply_gaussian_filter
        csq = imgaussfilt(csq, opts.gaussian_filter_sigma);
    end
elseif strcmpi(opts.version, 'pixi3')
    % read in image
    csq = double(imread(pfname));
    
    % ONLY CRRM MODE (PIXEL MODE) SUPPORTED - THIS IS FROM THE VENDOR
    % ct  = imread('/home/beams/S1IDUSER/mnt/s1a/misc/pixirad2/usb.after_repair/Calibrations/2010_crrm.tif');
    % ct  = fliplr(ct');
    ct  = load(opts.pfname_ct);
    ct  = ct.ct;
    
    bpt = load(opts.pfname_bpt);
    bpt = bpt.bpt;
    
    %%% SENSITIVITY CORRECTOIN
    if opts.apply_ct
        % csq = csq./ct;
        csq = csq.*ct;
    end
    
    %%% BAD PIXEL CORRECTION
    %%% MULTIPLE PASSES TO REMOVE ALL THE BAD PIXELS
    %%% MIGHT BE PATH DEPENDENT
    if opts.apply_bpt
        % [n, m]  = size(csq);
        
        there_are_bps   = true;
        counter = 0;
        while there_are_bps && (counter < 5)
            disp(sprintf('removing bad pixels - pass number %d', counter))
            disp(sprintf('removing bad pixels - %d bad pixels exist', sum(bpt(:))))
            [x_bpt, y_bpt]      = find(bpt == 1);
            
            %%% GET NEIGHBOR PIXEL VALUES
            for i = 1:1:length(x_bpt)
%                 pv  = nan(8,1);
%                 if (x_bpt(i)-1) >= 1
%                     if (y_bpt(i)-1) >= 1
%                         pv(1)   = csq(x_bpt(i)-1, y_bpt(i)-1); % nw
%                     end
%                     if (y_bpt(i)+1) <= m
%                         pv(2)   = csq(x_bpt(i)-1, y_bpt(i)+1); % ne
%                     end
%                     pv(3)   = csq(x_bpt(i)-1, y_bpt(i)); % n
%                 end
%                 
%                 if (x_bpt(i)+1) <= n
%                     if (y_bpt(i)-1) >= 1
%                         pv(4)   = csq(x_bpt(i)+1, y_bpt(i)-1); % sw
%                     end
%                     if (y_bpt(i)+1) <= m
%                         pv(5)   = csq(x_bpt(i)+1, y_bpt(i)+1); % se
%                     end
%                     pv(6)   = csq(x_bpt(i)+1, y_bpt(i)); % s
%                 end
%                 
%                 if (y_bpt(i)+1) <= m
%                     pv(7)   = csq(x_bpt(i), y_bpt(i)+1); % e
%                 end
%                 if (y_bpt(i)-1) >= 1
%                     pv(8)   = csq(x_bpt(i), y_bpt(i)-1); % w
%                 end
%                 pv(find(pv == 0))   = nan;
                % csq(x_bpt(i), y_bpt(i)) = mean(pv, 'omitnan');
                
                imc = imcrop(csq, [y_bpt(i)-1 x_bpt(i)-1 2 2]);
                imc = imc(:);
                imc(find(imc == 0)) = nan;
                
                %%% SOME USERS DO NOT HAVE OMITNAN
                if verLessThan('matlab', '9.0')
                    disp('this will work for now but upgrade matlab')
                    idx = ~isnan(imc);
                    csq(x_bpt(i), y_bpt(i)) = sum(imc(idx))/sum(idx);
                else
                    csq(x_bpt(i), y_bpt(i)) = mean(imc, 'omitnan');
                end
            end
            bpt = (csq == 0);
            
            there_are_bps = ~isempty(find(bpt == 1));
            
            counter = counter + 1;
        end
    end
end