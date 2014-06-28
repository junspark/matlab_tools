function pardata = ReadSpecParFile(fname, numchs)
% ReadSpecParFile - read par file generated from python
%
%   INPUT:
%
%   fname
%       name of the spec generated par file name
%
%   numchs
%       number of channels recorded in the spec generated par file (excluding the
%       date / time channel)
%
%   OUTPUT:
%
%   pardata
%       spec par file data in struct arrary.

fmtstring       = '%s %s %s %s %d %s';
for i = 1:1:numchs
    fmtstring   = [fmtstring, ' %f'];
end

fid     = fopen(fname, 'r','n');
textdata  = textscan(fid, fmtstring);
fclose(fid);

pardata.day     = textdata{1};
pardata.month   = textdata{2};
pardata.date    = textdata{3};
pardata.time    = textdata{4};
pardata.year    = textdata{5};
pardata.froot   = textdata{6};
pardata.chs     = zeros(length(pardata.date), numchs);
for i = 1:1:numchs
    pardata.chs(:,i)  = textdata{i+6};
end
