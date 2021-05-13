function pklist = GSAS2_Read_pkslst(pfname_pkslst)
% GSAS2_Read_pkslst - Reads GSAS2 pkslst file.
%
%   INPUT:
%
%   pfname_pkslst
%       full pathname to the pkslst file.
%
%   OUTPUT:
%
%   pklist
%       pklist parsed from pfname_pkslst file.

fid     = fopen(pfname_pkslst, 'r');

pklist  = [];
while ~feof(fid)
    linedata    = fgetl(fid);
    switch linedata(1)
        case '['
            pklist  = [pklist; eval(linedata)];
    end
end

fclose(fid);