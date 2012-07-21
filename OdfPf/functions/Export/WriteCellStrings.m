function WriteCellStrings(filename, strings)
% WRITECELLSTRINGS - Write a cell array of strings to a file.
%   
%   WriteCellStrings(filename, strings)
%
%   filename is a string, the name of the file
%   strings is a cell array of strings, the text to print
%
fid = fopen(filename, 'w');
fprintf(fid, '%s\n', strings{:});
fclose(fid);
