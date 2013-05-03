function imdata = ReadDXI11000(pfname, varargin)
% ReadDXI11000 - read DX-11000 file
%
%   INPUT:
%
%   pfname
%       name of the GE image stack file name
%
%   mode
%       data acquisition mode ('6fps' or '3fps')
%
%   OUTPUT:
%
%   imdata
%       image data in. for 6fps format the image is [2012 1340]. for 3fps
%       format the image is [4024 2680].
%       WHEN imagesc(imdata), TOP LEFT CORNER (0,0) IN THE FIGURE CORRESPONDS
%       TO THE TOP LEFT OF THE ACTIVE AREA ON THE SCINT-X DETECTOR WHEN THE
%       DETECTOR IS PLACED UP RIGHT.
%       CHECK W:\park_apr2013b\scintx\vff_scan_00019.tif AND
%       C:\Users\parkjs\Documents\Projects\VFF\park_apr2013b\setup\pics\DSC00291.JPG
%
%   NOTE:
%   1. if the image size changes, this file needs to be updated.

% default options
optcell = {...
    'Mode', '6fps', ...
    };

% update option
opts    = OptArgs(optcell, varargin);

imdata  = imread(pfname);

%%% IN THIS CONFIGURATION INCREASING RADIUS IS INCREASING COL NUMBER 
%%% INCREASING ETA IS INCREASING ROW NUMBER
imdata  = double(flipud(rot90(imdata, 2)));