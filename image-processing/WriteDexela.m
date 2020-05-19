function WriteDexela(filename, imdata)
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
%   imdata          image data to write out
%
%   OUTPUT:
%
%   none
disp(sprintf('writing \n %s', filename))
imwrite(imdata, filename, 'Compression', 'none');
