function Pixel = mm2Pixel(mm, c)
% mm2Pixel - converts mm to pixels
%
%   USAGE:
%
%   Pixel = mm2Pixel(mm)
%   Pixel = mm2Pixel(mm, c)
%
%   INPUT:
%
%   mm
%       length in mm
%
%   c
%       conversion factor between pixel and mm.  default is 100 micron = 1
%       pixel
%
%   OUTPUT:
%
%   pixel
%       length in pixel

if nargin == 1
    c   = 100*10^-6*1000;
end
Pixel  = mm/c;
