function BatchCorrection(path_bkg, bkg_num, root_bkg, ...
    path_image, root_image, ...
    genum, path_output, varargin)
% BatchCorrection - Corrects a batch of images. Currently implements only
% the background subtraction.
%
%   USAGE:
%
%   BatchCorrection(path_bkg, bkg_num, root_bkg, ...
%                       path_image, root_image, ...
%                       genum, path_output)
%   BatchCorrection(path_bkg, bkg_num, root_bkg, ...
%                       path_image, root_image, ...
%                       genum, path_output, varargin)
%
%   INPUT:
%
%   path_bkg            path where background file(s) are located.
%
%   bkg_num             if there are multiple background images, the low 
%                       and high values of the background images can be 
%                       provided and the function averages over all background images.
%
%   root_bkg            root of the background images
%
%	path_image          path where the the data image(s) are located.
%
%   root_image          root of the image(s) data
%	
%   genum               GE detector number (1 - 4)
%
%   path_output         output path
%
%   These arguments can be followed by a list of
%   parameter/value pairs which control certain plotting
%   features.  Options are:
%
%   CorrectAllImages    corrects all images in the path_image {false}. note
%                       that this option overrides any lo or hi image
%                       number inputs and root_image.
%
%   lo                  lowest number in the image series
%
%   hi                  highest number in the image series
%
%   OutAllFrames        output all corrected frames in a fastsweep file
%
%   NumDigits           number of digits in the file name (default = 6)
%
%   OUTPUT:             
%                   
%   (possibly) many files
%
%   NOTE:
%       1. Need to implement flipping (if needed)
%       2. This function is same as "BatchCorrection.py" (bad pixel
%       correction is slightly different)
%

% %%% START GCP
% delete(gcp); parpool(12);

%%% DEFAULT PARAMETERS
% npad        = '00000';
buffer_size = 8192;
frame_size  = 2048*2048*2;

% DEFAULT OPTIONS
optcell = {...
    'CorrectAllImages', true, ...
    'hi', -1, ...
    'lo', -1, ...
    'OutAllFrames', false, ...
    'DisplayFrames', false, ...
    'NumDigits', 6, ...
    };

% UPDATE OPTION
opts    = OptArgs(optcell, varargin);

%%% SETUP APPROPRIATE EXTENSION
ext_bkg     = ['ge', num2str(genum)];
ext_image   = ext_bkg;

%%% SETUP FILE NAME PATTERN
fname_fmt   = sprintf('%%s%%0%dd.%%s', opts.NumDigits);

%%% BAD PIXEL CORRECTION    %%% NOT YET IMPLEMENTED
BadPixelData    = LoadBadPixelData(genum);

%%% CREATE OUTPUT PATH
mkdir(path_output)

%%% CONSTRUCT BACKGROUND IMAGE
im_bkg  = zeros(2048,2048);
ct  = 0;
for i = 1:1:length(bkg_num)
    fname   = sprintf(fname_fmt, root_bkg, bkg_num(i), ext_bkg);
    % fname	= [root_bkg, npad(1:length(npad)-length(num2str(bkg_num(i)))) num2str(bkg_num(i)), '.', ext_bkg];
    pfname  = fullfile(path_bkg, fname);
    
    flist       = dir(pfname);
    if isempty(flist)
        disp(sprintf('check if %s exists', pfname));
        return
    end
    num_frame   = CalcNumFrames(flist.bytes, buffer_size, frame_size);
    
    for j = 1:1:num_frame
        im_bkg  = im_bkg + NreadGE(pfname, j);
        ct  = ct + 1;
    end
end
im_bkg  = im_bkg./ct;
if opts.DisplayFrames
    PlotImage(im_bkg, max(im_bkg(:)), min(im_bkg(:)))
end

if opts.CorrectAllImages
    cmd = ['flist = dir(''', path_image, filesep, '*.', ext_image, ''');'];
    eval(cmd)
    if isempty(flist)
        error('specified files do not exist');
    end
    
    parfor i = 1:length(flist)
        fname   = flist(i).name;
        pfname  = fullfile(path_image, fname);
        
        num_frame   = CalcNumFrames(flist(i).bytes, buffer_size, frame_size);
        disp(fname)
        
        sum_data    = zeros(2048,2048);
        for j = 1:1:num_frame
            disp('background subtraction in progress ...')
            frame_data  = NreadGE(pfname, j);
            if opts.OutAllFrames
                frame_data  = frame_data - im_bkg;
                frame_data  = CorrectBadPixels(frame_data, BadPixelData);
                
                fname_out   = [flist(i).name, '.frame.', num2str(j), '.cor'];
                pfname_out  = fullfile(path_output, fname_out);
                WriteSUM(pfname_out, frame_data);
            end
            sum_data    = sum_data + frame_data;
            
            if opts.DisplayFrames
                PlotImage(frame_data, max(frame_data(:)), min(frame_data(:)))
            end
        end
        if opts.DisplayFrames
            close all
        end
        
        %%% WRITE OUT SUM FILE
        if ~opts.OutAllFrames
            sum_data    = sum_data - num_frame.*im_bkg;
            sum_data    = CorrectBadPixels(sum_data, BadPixelData);
        end        
        fname_out   = [flist(i).name, '.sum'];
        pfname_out  = fullfile(path_output, fname_out);
        WriteSUM(pfname_out, sum_data);
        PlotImage(sum_data, max(sum_data(:)), min(sum_data(:)))
        title('Sum over all corrected frames')
        
        %%% WRITE OUT AVE FILE
        ave_data    = sum_data./num_frame;
        fname_out   = [flist(i).name, '.ave'];
        pfname_out  = fullfile(path_output, fname_out);
        WriteSUM(pfname_out, ave_data);
        PlotImage(ave_data, max(ave_data(:)), min(ave_data(:)))
        title('Average over all corrected frames')
    end
else
    image_num   = opts.lo:1:opts.hi;
    parfor i = 1:length(image_num)
        fname   = sprintf(fname_fmt, root_image, image_num(i), ext_bkg);
        % fname	= [root_image, npad(1:length(npad)-length(num2str(image_num(i)))) num2str(image_num(i)), '.', ext_image];
        pfname  = fullfile(path_bkg, fname);
        
        flist       = dir(pfname);
        num_frame   = CalcNumFrames(flist.bytes, buffer_size, frame_size);
        disp(fname);
        
        sum_data    = zeros(2048,2048);
        for j = 1:1:num_frame
            disp('background subtraction in progress ...')
            frame_data  = NreadGE(pfname, j);
            if opts.OutAllFrames
                frame_data  = frame_data - im_bkg;
                frame_data  = CorrectBadPixels(frame_data, BadPixelData);
                
                fname_out   = [flist.name, '.frame.', num2str(j), '.cor'];
                pfname_out  = fullfile(path_output, fname_out);
                WriteSUM(pfname_out, frame_data);
            end
            sum_data    = sum_data + frame_data;
            
            if opts.DisplayFrames
                PlotImage(frame_data, max(frame_data(:)), min(frame_data(:)))
            end
        end
        if opts.DisplayFrames
            close all
        end
        
        %%% WRITE OUT SUM FILE
        if ~opts.OutAllFrames
            sum_data    = sum_data - num_frame.*im_bkg;
            sum_data    = CorrectBadPixels(sum_data, BadPixelData);
        end
        fname_out   = [flist.name, '.sum'];
        pfname_out  = fullfile(path_output, fname_out);
        WriteSUM(pfname_out, sum_data);
        PlotImage(sum_data, max(sum_data(:)), min(sum_data(:)))
        title('Sum over all corrected frames')
        
        %%% WRITE OUT AVE FILE
        ave_data    = sum_data./num_frame;
        fname_out   = [flist.name, '.ave'];
        pfname_out  = fullfile(path_output, fname_out);
        WriteSUM(pfname_out, ave_data);
        PlotImage(ave_data, max(ave_data(:)), min(ave_data(:)))
        title('Average over all corrected frames')
    end
end
delete(gcp);