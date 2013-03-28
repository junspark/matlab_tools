function status = NwriteGE(pfname, data)
% NwriteGE - write GE file
%
%   INPUT:
%
%   pfname
%       name of the GE image full file name
%
%   data
%       data to write
%
%   OUTPUT:
%
%   img
%       image data in a [2048 x 2048] matrix
%
%   NOTE:
%   1. only writes 1 frame for now

fid     = fopen(pfname, 'w');
status  = fwrite(fid, zeros(4096,1), 'uint16');
status  = fwrite(fid, data, 'uint16');
status  = fclose(fid);