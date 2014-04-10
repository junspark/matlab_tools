function img = ReadSUM(filename)
% NreadGE - read GE file from fastsweep image stack
%
%   INPUT:
%
%   filename
%       name of the SUM or AVE image file name
%
%   OUTPUT:
%
%   img
%       image data in a [2048 x 2048] matrix

fp      = fopen(filename,'r','n');
img     = fread(fp,[2048 2048],'float');
fclose(fp);
