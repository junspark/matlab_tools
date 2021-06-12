function WriteDexela(filename, imdata, varargin)
% ReadDexela - reads .tiff images from Dexela.
%
%   USAGE:
%
%   WriteDexela(filename, imdata)
%
%   INPUT:
%
%   filename        File name of the Dexela tiff file.
%
%   imdata          image imdata to write out
%
%   OUTPUT:
%
%   none

disp(sprintf('writing \n %s', filename))
imwrite(imdata, filename, 'Compression', 'none');

% % update option
% optcell = {...
%     'apply_flip', true, ...
%     };
% 
% opts    = OptArgs(optcell, varargin);
% 
% if opts.apply_flip
%     imdata  = flipud(imdata);
% end
% 
% %%% MAKE IT INTO UINT16
% imdata  = uint16(imdata);
% 
% %%% START WRITING FILE
% tiff_class  = Tiff(filename, 'w');                                   %create object of Tiff class
% setTag(tiff_class, Tiff.TagID.ImageLength, size(imdata,1));          %define image dimentions
% setTag(tiff_class, Tiff.TagID.ImageWidth, size(imdata,2));
% 
% setTag(tiff_class, 'Photometric', Tiff.Photometric.MinIsBlack);    %define the color type of image
% 
% % specifies how image imdata components are stored on disk
% setTag(tiff_class, 'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
% 
% % Specify how to interpret each pixel sample (IEEEFP works with input doubles)
% setTag(tiff_class, 'BitsPerSample', 16);                              %because 1 double = 8byte = 64bits
% % setTag(tiff_class, 'SampleFormat', Tiff.SampleFormat.IEEEFP);
% 
% tiff_class.write(imdata);
% tiff_class.close;