function WriteDataFile(fname, text, data, basefmt, append)
% WriteDataFile - Write data to a file with a text header.
%   
%   USAGE:
%
%   WriteDataFile(fname, text, data)
%   WriteDataFile(fname, text, data, basefmt)
%   WriteDataFile(fname, text, data, basefmt, append)
%
%   INPUT:
%
%   fname is a string, 
%         the name of the file
%   text  is a cell array of strings, 
%         the text header to write at the beginning of the file; 
%         it can be empty
%   data  is m x n, (numeric)
%         the data to write; the data is written as passed, 
%         m rows of n columns each
%   basefmt is a string, (optional, default = '%25.16e')
%         the printf-style format to use on each value;  if this 
%         argument is left out or is the empty string '', the default
%         value for reals is used; if the argument is
%         the string 'integer', the format used is '%10d';
%         other values are used without change
%   append is a scalar,  (optional, default = 0)
%          if nonzero, the output file is appended to instead of 
%          overwritten
%
%   OUTPUT:  none
%
[m, n] = size(data);
%
dflt_fmt = '%25.16e';
dflt_int = '%10d';
%
if (nargin < 4)
  basefmt = dflt_fmt;
end
%
if (strcmp(basefmt, 'integer'))
  basefmt = dflt_int;
elseif (strcmp(basefmt, ''))
  basefmt = dflt_fmt;
end
%
if (nargin < 5)
  append = 0;
end
%
if (append)
  openfor = 'a';
else
  openfor = 'w';
end
%
%  Construct format for printf.
%
fmtstr = basefmt;
for i=2:n
  fmtstr = [fmtstr, ' ', basefmt];
end
fmtstr = [fmtstr, '\n'];
%
%  Open file and write text and data.
%
f = fopen(fname, openfor); 
if (f == -1)
  error('failed to open file')
end
%
fprintf(f, '%s\n', text{:});
fprintf(f, fmtstr, data');  % so that it is written in correct order
%
fclose(f);
%
