function ExportDX_PointData(fname, pos, data)
% ExportDX_PointData - export point data into DX format
%
%   USAGE:
%
%   ExportDX_PointData(fname, pos, data)
%
%   INPUT:
%
%   fname
%       output dx-file name
%
%   pos
%       coordinates of the points (n x m)
%
%   data
%       corresponding data for the points (n x 1)
%
%   OUTPUT:
%
%   none

fname_pos   = [fname, '.posdat'];

[np, mp]    = size(pos);
[nd, md]    = size(data);

if np ~= nd
    error('number of points does not match the number of data')
    return
end

data    = [pos; data];
ndata   = size(data,1);

PrnStr  = [];
for i = 1:1:ndata
    PrnStr  = [PrnStr, '%f '];
end
PrnStr  = [PrnStr, '\n'];

fid_pos = fopen(fname_pos, 'w');
fprintf(fid_pos, PrnStr, data);
fclose(fid_pos);