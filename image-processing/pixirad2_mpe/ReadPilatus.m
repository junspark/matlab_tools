function [csq] = ReadPilatus(pfname, varargin)
% ReadPilatus - read Pilatus pattern file.
%
%   INPUT:
%
%   pfname
%       name of the Pixirad image file.
%
%   format
%       tiff or raw format output. Default = tiff.
%
%   OUTPUT:
%
%   csq
%       Pixirad image file information mapped to an image with square
%       pixels

% default options
optcell = {...
    'format', 'tiff', ...
    };
opts    = OptArgs(optcell, varargin);

if strcmpi(opts.format, 'tiff')
    disp('reading tiff pilatus file.')
    csq = double(imread(pfname));
elseif strcmpi(opts.format, 'tiff')
    error('raw pilatus reader not implemented.')
else
    error('unknown format.')
end