function ExportDX_MeshData(fname, pos, data, con)
% ExportDX_MeshData - export mesh data into DX format
%
%   INPUT:
%
%   fname
%       output dx-file name
%
%   pos
%       coordinates of the mesh points
%
%   data
%       corresponding data for the points
%
%   con
%       connectivity of the mesh
%
%   OUTPUT:
%
%   none

fname_pos   = [fname, '.pos'];
fname_con   = [fname, '.con'];

fid_pos = fopen(fname_pos, 'w');
fid_con = fopen(fname_con, 'w');

data    = [pos; data];
ndata   = size(data,1);

PrnStr  = [];
for i = 1:1:ndata
    PrnStr  = [PrnStr, '%f '];
end
PrnStr  = [PrnStr, '\n'];
fprintf(fid_pos, PrnStr, data);

ncon    = size(con, 1);
PrnStr  = [];
for i = 1:1:ncon
    PrnStr  = [PrnStr, '%d '];
end
PrnStr  = [PrnStr, '\n'];
fprintf(fid_con, PrnStr, con);

fclose(fid_pos);
fclose(fid_con);