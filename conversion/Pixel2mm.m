function mm = Pixel2mm(Pixel, c)
% Pixel2mm - converts pixels to mm
%
%   USAGE:
%
%   mm = Pixel2mm(Pixel)
%   mm = Pixel2mm(Pixel, c)
%
%   INPUT:
%
%   Pixel
%       length in Pixel
%
%   c
%       conversion factor between pixel and mm.  default is 100 micron = 1
%       pixel
%
%   OUTPUT:
%
%   mm
%       length in mm

if nargin == 1
    c   = 100*10^-6*1000;
end
mm  = Pixel*c;