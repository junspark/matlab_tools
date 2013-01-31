function FileData = ReadEnduratecData(FileName, NumCols)
% ReadEnduratecData - import text file exported from Enduratec load frame
%
%   USAGE:
%
%   FileData = ReadEnduratecData(FileName, NumCols)
%
%   INPUT:
%
%   FileName
%       name of the Enduractec text file (not its proprietary binary files)
%
%   NumCols
%       number of columns (fields) in the text file
%
%   OUTPUT:
%
%   FileData
%       tabular data from the Enduractec text file
%

%%%
%clear all
%close all
%clc
%FileName    = 'C:\mpm_research\MaterialProperty\Ti64\Navy.Ti64\Ti.Navy.TOut.svse.txt';

fid     = fopen(FileName, 'r');
FullFileName    = fgetl(fid);
fseek(fid, 2, 'cof');

FileData.DataLabel  = fgetl(fid);
FileData.DataUnits  = fgetl(fid);
FileData.Data       = fscanf(fid, '%f,');

NumRows = length(FileData.Data);
NumRows = NumRows/NumCols;

FileData.Data   = reshape(FileData.Data, NumCols, NumRows)';
fclose(fid);