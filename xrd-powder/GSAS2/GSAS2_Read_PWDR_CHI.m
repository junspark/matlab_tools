function data_chi = GSAS2_Read_PWDR_CHI(pfname_chi)
% GSAS2_Read_PWDR_CHI - GSAS2 PWDR CHI data
%
%   INPUT:
%
%   pfname_chi
%       full pathname of the GSAS2 PWDR CHI data file.
%
%   OUTPUT:
%
%   data_chi
%       tabular data of the the GSAS2 PWDR CHI data file.


data_chi    = table2array(readtable(pfname_chi, 'FileType', 'text', ...
    'ReadVariableNames', false, ...
    'HeaderLines', 4));