function img = NreadDexela(pfname, frameno)
% NreadDexela - read Dexela file from fastsweep image stack saved in ncdf
% format (nc extension).
%
%   INPUT:
%
%   filename
%       name of the Dexela image stack file name
%
%   frameno
%       frame to load
%
%   OUTPUT:
%
%   img
%       image data in a [3072 x 3888] matrix in default size - 1 x 1 binning
%

% finfo       = ncinfo(pfname);
% nframe      = finfo.Dimensions(1).Length;
% numpixHorz  = finfo.Dimensions(2).Length;
% numpixVert  = finfo.Dimensions(3).Length;
% 
% if frameno > nframe
%     error('frame number %d requested is larger than number of frames %d stored in the file %s', frameno, nframe, pfname);
% else
%     istart  = [1 1 frameno]';
%     iend    = [numpixVert numpixHorz 1]';
%     img     = ncread(pfname, 'array_data', istart, iend);
% end
% 
% img = double(img);

%%% WE SAVE IN TIF FOR DEXELA
img = imread(pfname);
img = double(img);