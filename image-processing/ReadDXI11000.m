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
imdata  = double(rot90(imdata, 2)); 