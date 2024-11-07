function BatchCorrection_varex(path_bkg, bkg_num, root_bkg, ...
    path_image, root_image, ...
    path_output, varargin)
% BatchCorrection_varex - Corrects a batch of image for Varex. HDF5 files
% only.
%
%   USAGE:
%
%   BatchCorrection_varex(path_bkg, bkg_num, root_bkg, ...
%                   path_image, root_image, path_output)
%
%   BatchCorrection_varex(path_bkg, bkg_num, root_bkg, ...
%                   path_image, root_image, genum, path_output, varargin)
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
%   FramesToIgnore      Frames to ignore (to ignore shadowed frames) 
%                       (default = nan). Otherwise, this is a {n x 1} cell 
%                       array where each cell element contains a [m x 1] 
%                       matrix indicating which frames to ignore.
%
%   SumOnly             Only output sum files. If 0, ave file will be saved
%                       as well (default: 1).
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

% DEFAULT OPTIONS
optcell = {...
    'CorrectAllImages', true, ...
    'hi', -1, ...
    'lo', -1, ...
    'OutAllFrames', false, ...
    'DisplayFrames', false, ...
    'NumDigits', 6, ...
    'FramesToIgnore', 'none', ...
    'SumOnly', 1, ...
    };

% UPDATE OPTION
opts    = OptArgs(optcell, varargin);

%%% SETUP APPROPRIATE EXTENSION
ext_bkg     = 'vrx.h5';
ext_image   = ext_bkg;

fname_fmt   = sprintf('%%s_%%0%dd.%%s', opts.NumDigits);

%%% BAD PIXEL CORRECTION
% BadPixelData    = LoadBadPixelData(genum);

%%% CREATE OUTPUT PATH
if isfolder(path_output) == 0
    mkdir(path_output)
end

%%% CONSTRUCT BACKGROUND IMAGE
im_bkg  = zeros(2880,2880);
ct  = 0;
for iii = 1:1:length(bkg_num)
    fname   = sprintf(fname_fmt, root_bkg, bkg_num(iii), ext_bkg);
    pfname  = fullfile(path_bkg, fname);
    
    flist       = dir(pfname);
    if isempty(flist)
        error('check if %s exists\n', pfname);
    end
    
    if iii == 1
        im_bkg      = zeros(h5info(pfname, '/exchange/bright/').ChunkSize(1), h5info(pfname, '/exchange/bright/').ChunkSize(2));
    end
    num_frames  = h5info(pfname, '/exchange/bright/').Dataspace.Size(3);
    
    im_bkg_iii  = h5read(pfname, '/exchange/bright/');
    im_bkg_iii  = sum(im_bkg_iii,3);
    im_bkg      = im_bkg + im_bkg_iii;
    ct          = ct + num_frames;
end
im_bkg  = im_bkg./ct;

% figure,
% imagesc(im_bkg_iii)
% 
% figure,
% imagesc(im_bkg)
% return
%%% DISPLAY BACKGROUND IF REQUESTED
if opts.DisplayFrames
    PlotImage(im_bkg, max(im_bkg(:)), min(im_bkg(:)))
end

%%% DETERMINE FRAMES TO CORRECT
if strcmpi(opts.FramesToIgnore, 'none')
    disp('all frames in the image stack will be corrected ...');
    FramesToIgnore  = nan;
else
    return
    disp('some frames will be ignored ...');
    FramesToIgnore  = [opts.FramesToIgnore{:}];
    FramesToIgnore  = unique(FramesToIgnore);
    idx = FramesToIgnore <= 0;
    if sum(idx) > 0
        %disp('check frame numbers to ingore ...');
        %disp('positive integers only ...');
        %return
        error('check frame numbers to ignore ...\npositive integers only ...\n'); % edited by Connor Horn 7/31/17
    end
end

%%% DETERMINE FILES TO CORRECT
if opts.CorrectAllImages
    %%% CORRECT ALL IMAGES
    %%% GET A LIST OF ALL IMAGES MATCHING THE PATTERN
    cmd = ['flist = dir(''', path_image, filesep, root_image, '*.', ext_image, ''');'];
    
    eval(cmd)
else
    lo_image_num    = opts.lo;
    hi_image_num    = opts.hi;
    
    if hi_image_num < lo_image_num
        error('hi needs to be larger than lo');
    end
    
    image_num       = lo_image_num:1:hi_image_num;
    
    ct  = 1;
    for iii = 1:1:length(image_num)
        fname   = sprintf(fname_fmt, root_image, image_num(iii), ext_image);
        pfname  = fullfile(path_image, fname);
        
        if isfile(pfname)
            flist(ct,1) = dir(pfname);
            ct  = ct + 1;
        end
    end
end

if isempty(flist)
    error('specified files do not exist');
end

parfor iii = 1:1:length(flist)
    fname   = flist(iii).name;
    pfname  = fullfile(path_image, fname);
    
    sum_data    = h5read(pfname, '/exchange/bright/');
    sum_data    = mean(sum_data,3);
    sum_data    = sum_data - im_bkg;
    
    if opts.DisplayFrames
        PlotImage(sum_data, max(sum_data(:)), min(sum_data(:)))
    end
    
    fname_out   = [flist(iii).name, '.sum.tiff'];
    pfname_out  = fullfile(path_output, fname_out);
    WriteVarexTIFF(pfname_out, sum_data)
    
    fname_out   = [flist(iii).name, '.sum.h5'];
    pfname_out  = fullfile(path_output, fname_out);
    % WriteVarexHDF5(filename, imdata, metadata, varargin)
end

%     flist
%     return   %%% NEED TO CHECK BELOW
%     
%     %%% GO THROUGH THE FILE LIST
%     parfor i = 1:length(flist)
%         fname   = flist(i).name;
%         pfname  = fullfile(path_image, fname);
%         
%         %%% DETERMINE NUMBER OF FRAMES IN THIS STACK
%         num_frames   = CalcNumFrames(flist(i).bytes, buffer_size, frame_size);
%         
%         %%% CHECK IF THE FRAMES TO SKIP MAKES SENSE
%         if max(FramesToIgnore) > num_frames
%             %%% SKIP IF MAX FRAME NUMBER TO IGNORE IS LARGER THAN THE TOTAL
%             %%% NUMBER OF FRAMES IN THEH STACK
%             disp('OH NO! The number of frames in the image stack inconsistent with frames requested for correction ...')
%             fprintf('skipping %s ...\n', fname);
%         else
%             %%% CONTINUE IF THE MAX FRAME NUMBER TO IGNORE IS SMALLER THAN THE TOTAL
%             %%% NUMBER OF FRAMES IN THE IMAGE STACK
%             FramesToCorrect     = 1:1:num_frames;
%             if ~strcmpi(opts.FramesToIgnore, 'none')
%                 for j = 1:1:length(FramesToIgnore)
%                     idx = FramesToCorrect == FramesToIgnore(j);
%                     FramesToCorrect     = FramesToCorrect(~idx);
%                 end
%             end
%             
%             sum_data    = zeros(2048,2048);
%             for j = FramesToCorrect
%                 fprintf('%d frame of %s - background subtraction in progress ...\n', j, fname);
%                 frame_data  = NreadGE(pfname, j);
%                 if opts.OutAllFrames
%                     frame_data  = frame_data - im_bkg;
%                     frame_data  = CorrectBadPixels(frame_data, BadPixelData);
% 
%                     % fname_out   = [flist(i).name, '.frame.', num2str(j), '.cor'];
%                     fname_out   = [flist(i).name, '_frame_', num2str(j), '.cor32'];
%                     pfname_out  = fullfile(path_output, fname_out);
%                     WriteSUM(pfname_out, frame_data);
%                 end
%                 sum_data    = sum_data + frame_data;
% 
%                 if opts.DisplayFrames
%                     PlotImage(frame_data, max(frame_data(:)), min(frame_data(:)))
%                 end
%             end
%             if opts.DisplayFrames
%                 close all
%             end
% 
%             %%% WRITE OUT SUM FILE
%             if ~opts.OutAllFrames
%                 sum_data    = sum_data - im_bkg*length(FramesToCorrect);
%                 sum_data    = CorrectBadPixels(sum_data, BadPixelData);
%             end
%             fname_out   = [flist(i).name, '.sum'];
%             pfname_out  = fullfile(path_output, fname_out);
%             WriteSUM(pfname_out, sum_data);
%             
%             if opts.DisplayFrames
%                 PlotImage(sum_data, max(sum_data(:)), min(sum_data(:)))
%                 title('Sum over all corrected frames')
%             end
% 
%             %%% WRITE OUT AVE FILE
%             ave_data    = sum_data./length(FramesToCorrect);
%             fname_out   = [flist(i).name, '.ave'];
%             pfname_out  = fullfile(path_output, fname_out);
%             if ~opts.SumOnly
%                 WriteSUM(pfname_out, ave_data);
%             end
% 
%             if opts.DisplayFrames
%                 PlotImage(ave_data, max(ave_data(:)), min(ave_data(:)))
%                 title('Average over all corrected frames')
%             end
%         end
%     end
% else
%     %%% SELECTIVELY CORRECT IMAGES
%     lo_image_num = opts.lo;
%     hi_image_num = opts.hi;
%     image_num   = lo_image_num:1:hi_image_num;
%     
%     if image_num == -1
%         error('correction did not run\need to specify image numbers\n'); % edited by CH 7/31/17
%     end
%     
%     parfor i = 1:length(image_num)
%         fname   = sprintf(fname_fmt, root_image, image_num(i), ext_image);
%         pfname  = fullfile(path_bkg, fname);
%         
%         if isfile(pfname)
%             fprintf('working on %s\n',pfname)
%             flist       = dir(pfname);
%             
%             %%% DETERMINE NUMBER OF FRAMES IN THIS STACK
%             num_frames   = CalcNumFrames(flist(1).bytes, buffer_size, frame_size);
%             
%             %%% CHECK IF THE FRAMES TO SKIP MAKES SENSE
%             if max(FramesToIgnore) > num_frames
%                 %%% SKIP IF MAX FRAME NUMBER TO IGNORE IS LARGER THAN THE TOTAL
%                 %%% NUMBER OF FRAMES IN THEH STACK
%                 disp('OH NO! The number of frames in the image stack inconsistent with frames requested for correction ...')
%                 fprintf('skipping %s ...\n', fname);
%             else
%                 %%% CONTINUE IF MAX FRAME NUMBER TO IGNORE IS SMALLER THAN THE TOTAL
%                 %%% NUMBER OF FRAMES IN THEH STACK
%                 FramesToCorrect     = 1:1:num_frames;
%                 if ~strcmpi(opts.FramesToIgnore, 'none')
%                     for j = 1:1:length(FramesToIgnore)
%                         idx = FramesToCorrect == FramesToIgnore(j);
%                         FramesToCorrect     = FramesToCorrect(~idx);
%                     end
%                 end
%                 
%                 sum_data    = zeros(2048,2048);
%                 for j = FramesToCorrect
%                     fprintf('%d frame of %s - background subtraction in progress ...\n', j, fname);
%                     frame_data  = NreadGE(pfname, j);
%                     if opts.OutAllFrames
%                         frame_data  = frame_data - im_bkg;
%                         frame_data  = CorrectBadPixels(frame_data, BadPixelData);
%                         
%                         % fname_out   = [flist.name, '.frame.', num2str(j), '.cor'];
%                         fname_out   = [flist.name, '_frame_', num2str(j), '.cor32'];
%                         pfname_out  = fullfile(path_output, fname_out);
%                         WriteSUM(pfname_out, frame_data);
%                     end
%                     sum_data    = sum_data + frame_data;
%                     
%                     if opts.DisplayFrames
%                         PlotImage(frame_data, max(frame_data(:)), min(frame_data(:)))
%                     end
%                 end
%                 if opts.DisplayFrames
%                     close all
%                 end
%                 
%                 %%% WRITE OUT SUM FILE
%                 if ~opts.OutAllFrames
%                     sum_data    = sum_data - im_bkg*length(FramesToCorrect);
%                     sum_data    = CorrectBadPixels(sum_data, BadPixelData);
%                 end
%                 fname_out   = [flist.name, '.sum'];
%                 pfname_out  = fullfile(path_output, fname_out);
%                 WriteSUM(pfname_out, sum_data);
%                 
%                 if opts.DisplayFrames
%                     PlotImage(sum_data, max(sum_data(:)), min(sum_data(:)))
%                     title('Sum over all corrected frames')
%                 end
%                 
%                 %%% WRITE OUT AVE FILE
%                 ave_data    = sum_data./length(FramesToCorrect);
%                 fname_out   = [flist.name, '.ave'];
%                 pfname_out  = fullfile(path_output, fname_out);
%                 if ~opts.SumOnly
%                     WriteSUM(pfname_out, ave_data);
%                 end
%                 
%                 if opts.DisplayFrames
%                     PlotImage(ave_data, max(ave_data(:)), min(ave_data(:)))
%                     title('Average over all corrected frames')
%                 end
%             end
%         else
%             fprintf('%s does not exist. skipping.\n',pfname)
%         end
%     end
% end
% 
% %%% END GCP
% if license('test', 'distrib_computing_toolbox') && isunix
%     disp(sprintf('parallel computing toolbox available'));
%     delete(gcp);
% else
%     disp(sprintf('parallel computing toolbox unavailable'));
% end
