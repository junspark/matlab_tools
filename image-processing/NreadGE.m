function img = NreadGE(filename, frameno)
% NreadGE - read GE file from fastsweep image stack
%
%   INPUT:
%
%   filename
%       name of the GE image stack file name
%
%   frameno
%       frame to load
%
%   OUTPUT:
%
%   img
%       image data in a [2048 x 2048] matrix
%
%   NOTE:
%   1. if the image size changes, this file needs to be updated.
%   2. if the buffer contains useful information, this file needs to read
%   it in.

buffer  = 8192;

fp      = fopen(filename,'r','n');
offset  = buffer + (frameno-1)*2048*2048*2;

fseek(fp,offset,'bof');

img     = fread(fp,[2048 2048],'*uint16');

fclose(fp);

img = double(img);
