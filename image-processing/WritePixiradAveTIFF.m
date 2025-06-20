function [] = WritePixiradAveTIFF(pfname, data, varargin)
% WritePixiradAveTIFF(pfname, data) - write pixirad ave uint16 tiff 
% readable by GSAS2
%
%   INPUT:
%
%   pfname
%       name of the image full file name
%
%   data
%       data to write
%
%   OPTIONAL INPUT:
%
%   apply_flip (default: true)
%       apply flip about horizontal plane. 
% 
%   bit_depth (default: single)
%       output pixel can be "single" or "uint16".
%
%   OUTPUT:
%
%   status
%       status of the file write
%

% update option
optcell = {...
    'apply_flip', true, ...
    'bit_depth', 'single', ...
    };

opts    = OptArgs(optcell, varargin);

if opts.apply_flip
    data    = flipud(data);
end

if strcmp(opts.bit_depth, 'uint16')
    disp(sprintf('*** writing %s as %s', pfname, opts.bit_depth))

    %%% MAKE IT INTO UINT16
    data   = uint16(data);

    %%% START WRITING FILE
    tiff_class  = Tiff(pfname, 'w');                                   %create object of Tiff class
    setTag(tiff_class, Tiff.TagID.ImageLength, size(data,1));          %define image dimentions
    setTag(tiff_class, Tiff.TagID.ImageWidth, size(data,2));

    setTag(tiff_class, 'Photometric', Tiff.Photometric.MinIsBlack);    %define the color type of image

    % specifies how image data components are stored on disk
    setTag(tiff_class, 'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);

    % Specify how to interpret each pixel sample (IEEEFP works with input doubles)
    setTag(tiff_class, 'BitsPerSample', 16);                              %because 1 double = 8byte = 64bits
    % setTag(tiff_class, 'SampleFormat', Tiff.SampleFormat.IEEEFP);

    tiff_class.write(data);
    tiff_class.close;
    
elseif strcmp(opts.bit_depth, 'single')
    disp(sprintf('*** writing %s as %s', pfname, opts.bit_depth))

    %%% MAKE IT INTO SINGLE
    data   = single(data);

    tiff_class  = Tiff(pfname, 'w');                                   %create object of Tiff class
    setTag(tiff_class, Tiff.TagID.ImageLength, ...
        size(imdata_out_per_sec,1));          %define image dimentions
    setTag(tiff_class, Tiff.TagID.ImageWidth, ...
        size(imdata_out_per_sec,2));

    setTag(tiff_class, 'Photometric', Tiff.Photometric.MinIsBlack);    %define the color type of image

    % specifies how image data components are stored on disk
    setTag(tiff_class, 'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);

    % Specify how to interpret each pixel sample (IEEEFP works with input doubles)
    setTag(tiff_class, 'BitsPerSample', 32);                              %because 1 double = 8byte = 64bits
    setTag(tiff_class, 'SampleFormat', Tiff.SampleFormat.IEEEFP);

    tiff_class.write(data);
    tiff_class.close;
else
    warning(sprintf('*** format %s does not exist', opts.bit_depth))
end