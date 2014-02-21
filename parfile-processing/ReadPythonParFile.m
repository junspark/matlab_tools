function pardata = ReadPythonParFile(fname, numchs)
% ReadPythonParFile - read par file generated from python
%
%   INPUT:
%
%   fname
%       name of the python generated par file name
%
%   numchs
%       number of channels recorded in the python generated par file (excluding the
%       date / time channel)
%
%   OUTPUT:
%
%   pardata
%       strip chart data in struct arrary.
fmtstring       = '%s %s %d %s %d';
for i = 1:1:numchs
    fmtstring   = [fmtstring, ', %f'];
end

fid     = fopen(fname, 'r','n');
header  = fgetl(fid);
textdata  = textscan(fid, fmtstring);
fclose(fid);

pardata.hdr     = header;
pardata.day     = textdata{1};
pardata.month   = textdata{2};
pardata.date    = textdata{3};
pardata.time    = textdata{4};
pardata.year    = textdata{5};
pardata.chs     = zeros(length(pardata.date), numchs);
for i = 1:1:numchs
    pardata.chs(:,i)  = textdata{i+5};
end