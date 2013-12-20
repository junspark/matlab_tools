function scdata = ReadStripChartData(filename, numchs)
% ReadStripChartData - read strip chart data file
%
%   INPUT:
%
%   filename
%       name of the strip chart file name
%
%   numchs
%       number of channels recorded in the strip chart file (excluding the
%       date / time channel)
%
%   OUTPUT:
%
%   scdata
%       strip chart data in struct arrary.
fmtstring       = '%s %s';
for i = 1:1:numchs
    fmtstring   = [fmtstring, ' %s'];
end

fid     = fopen(filename, 'r','n');
header  = fgetl(fid);
textdata  = textscan(fid, fmtstring);
fclose(fid);

idxBadVal   = find(strcmp(textdata{end}, 'BadVal'));
if ~isempty(idxBadVal)
    for i = 1:1:(numchs+2)
        textdata{i}   = textdata{i}(idxBadVal(end):end);
    end
end

scdata.hdr  = header;
scdata.date = textdata{1};
scdata.time = textdata{2};
scdata.chs  = zeros(length(scdata.date), numchs);
for i = 1:1:numchs
    scdata.chs(:,i)  = str2double(textdata{i+2});
end