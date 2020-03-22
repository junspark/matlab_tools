function data = Read_s6bmb_edf(filename)
% Read_s6bmb_edf - reads .edf for S6BMB calibration information.
%
%   USAGE:
%
%   data = Read_s6bmb_edf(filename)
%
%   INPUT:
%
%   filename    File name of the Dexela tiff file.
%
%   OUTPUT:
%
%   data        struct with polynominal and tth information for
%               each detector from detector 1 - 10.

data_tbl    = readtable(filename, ...
    'FileType', 'text', 'ReadVariableNames', 0);
data_tbl    = reshape(data_tbl.Variables, 6, 10);

%%% DET ORDER detnum  = [1 2 3 4 5 6 7 8 9 0];
for i = 1:1:10
    for j = 1:1:4
        data(i).polynomial(j)   = str2double(data_tbl{j,i});
    end
    data(i).tth = str2double(data_tbl{5,i});
end