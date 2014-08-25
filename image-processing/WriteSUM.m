function status = WriteSUM(pfname, data)
% WriteSUM - write SUM file
%
%   INPUT:
%
%   pfname
%       name of the SUM image full file name
%
%   data
%       data to write. 2048 x 2048
%
%   OUTPUT:
%
%   status
%       status of the file write
%

fid = fopen(pfname, 'w');
status  = fwrite(fid, data, 'float');
status  = fclose(fid);