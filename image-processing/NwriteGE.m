function status = NwriteGE(pfname, data)
% NwriteGE - write GE file
%
%   INPUT:
%
%   pfname
%       name of the GE image full file name
%
%   data
%       data to write. 2048 x 208 x n where n is the number of frames
%
%   OUTPUT:
%
%   status
%       status of the file write
%

nFrames = size(data, 3);
fid = fopen(pfname, 'w');
status  = fwrite(fid, zeros(4096,1), 'uint16');
for i = 1:1:nFrames
    status  = fwrite(fid, data(:,:,i), 'uint16');
end
status  = fclose(fid);