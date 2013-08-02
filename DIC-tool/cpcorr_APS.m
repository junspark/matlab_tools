function xyinput = cpcorr_APS(varargin)
%CPCORR Tune control point locations using cross-correlation. 
%   INPUT_POINTS = CPCORR(INPUT_POINTS_IN,BASE_POINTS_IN,INPUT,BASE) uses
%   normalized cross-correlation to adjust each pair of control points
%   specified in INPUT_POINTS_IN and BASE_POINTS_IN.
%
%   INPUT_POINTS_IN must be an M-by-2 double matrix containing the
%   coordinates of control points in the input image.  BASE_POINTS_IN is
%   an M-by-2 double matrix containing the coordinates of control points
%   in the base image.
%
%   CPCORR returns the adjusted control points in INPUT_POINTS, a double
%   matrix the same size as INPUT_POINTS_IN.  If CPCORR cannot correlate a
%   pairs of control points, INPUT_POINTS will contain the same coordinates
%   as INPUT_POINTS_IN for that pair.
%
%   CPCORR will only move the position of a control point by up to 4
%   pixels.  Adjusted coordinates are accurate up to one tenth of a
%   pixel.  CPCORR is designed to get subpixel accuracy from the image
%   content and coarse control point selection.
%
%   Note that the INPUT and BASE images must have the same scale for
%   CPCORR to be effective.
%
%   CPCORR cannot adjust a point if any of the following occur:
%     - points are too near the edge of either image
%     - regions of images around points contain Inf or NaN
%     - region around a point in input image has zero standard deviation
%     - regions of images around points are poorly correlated
%
%   Class Support
%   -------------
%   The images can be numeric and must contain finite values. The input
%   control point pairs are double.
%
%   Example
%   --------
%   This example uses CPCORR to fine-tune control points selected in an
%   image.  Note the difference in the values of the INPUT_POINTS matrix
%   and the INPUT_POINTS_ADJ matrix.
%
%       input = imread('onion.png');
%       base = imread('peppers.png');
%       input_points = [127 93; 74 59];
%       base_points = [323 195; 269 161];
%       input_points_adj = cpcorr(input_points,base_points,...
%                                 input(:,:,1),base(:,:,1))
%
%   See also CP2TFORM, CPSELECT, NORMXCORR2, IMTRANSFORM.
%
%   NOTE: THIS FUNCTION WAS REVISED FOR STRAIN CALCULATION APPLICATION; THE
%   ORIGINAL MATLAB FUNCTION IS "cpcorr".

%   Copyright 1993-2011 The MathWorks, Inc.
%   $Revision: 1.16.4.8.2.1 $  $Date: 2011/07/18 00:33:04 $

%   Input-output specs
%   ------------------
%   INPUT_POINTS_IN: M-by-2 double matrix 
%              INPUT_POINTS_IN(:)>=0.5
%              INPUT_POINTS_IN(:,1)<=size(INPUT,2)+0.5
%              INPUT_POINTS_IN(:,2)<=size(INPUT,1)+0.5
%
%   BASE_POINTS_IN: M-by-2 double matrix 
%              BASE_POINTS_IN(:)>=0.5
%              BASE_POINTS_IN(:,1)<=size(BASE,2)+0.5
%              BASE_POINTS_IN(:,2)<=size(BASE,1)+0.5
%
%   INPUT:   2-D, real, full matrix
%            logical, uint8, uint16, or double
%            must be finite (no NaNs, no Infs inside regions being correlated)
%
%   BASE:    2-D, real, full matrix
%            logical, uint8, uint16, or double
%            must be finite (no NaNs, no Infs inside regions being correlated)

[xyinput_in,xybase_in,input,base] = ParseInputs(varargin{:});

CORRSIZE = 15;

% get all rectangle coordinates
rects_input = calc_rects(xyinput_in,CORRSIZE,input);
rects_base = calc_rects(xybase_in,2*CORRSIZE,base);

ncp = size(xyinput_in,1);

xyinput = xyinput_in; % initialize adjusted control points matrix

for icp = 1:ncp

    if isequal(rects_input(icp,3:4),[0 0]) || ...
       isequal(rects_base(icp,3:4),[0 0]) 
        % near edge, unable to adjust
        continue
    end
    
    sub_input = imcrop(input,rects_input(icp,:));
    sub_base = imcrop(base,rects_base(icp,:));    

    inputsize = size(sub_input);

    % make sure finite
    if any(~isfinite(sub_input(:))) || any(~isfinite(sub_base(:)))
        % NaN or Inf, unable to adjust
        continue
    end

    % check that template rectangle sub_input has nonzero std
    if std(sub_input(:))==0
        % zero standard deviation of template image, unable to adjust
        continue
    end

    norm_cross_corr = normxcorr2(sub_input,sub_base);    

    % get subpixel resolution from cross correlation
    subpixel = true;
    [xpeak, ypeak, amplitude] = findpeak_APS(norm_cross_corr,subpixel);

    % eliminate any poor correlations
    THRESHOLD = 0.5;
    if (amplitude < THRESHOLD) 
        % low correlation, unable to adjust
        continue
    end
    
    % offset found by cross correlation
    corr_offset = [ (xpeak-inputsize(2)-CORRSIZE) (ypeak-inputsize(1)-CORRSIZE) ];

    % eliminate any big changes in control points
    ind = find(abs(corr_offset) > (CORRSIZE-1), 1);
    if ~isempty(ind)
        % peak of norxcorr2 not well constrained, unable to adjust
        continue
    end

    % Change for digital image correlation code
    % input_fractional_offset = xyinput(icp,:) - round(xyinput(icp,:));
    % base_fractional_offset = xybase_in(icp,:) - round(xybase_in(icp,:));    

    input_fractional_offset = xyinput(icp,:) - round(xyinput(icp,:)*1000)/1000;
    base_fractional_offset = xybase_in(icp,:) - round(xybase_in(icp,:)*1000)/1000; 
    
    % adjust control point
    xyinput(icp,:) = xyinput(icp,:) - input_fractional_offset - corr_offset + base_fractional_offset;

end

%-------------------------------
%
function rect = calc_rects(xy,halfwidth,img)

% Calculate rectangles so imcrop will return image with xy coordinate inside center pixel

default_width = 2*halfwidth;
default_height = default_width;

% xy specifies center of rectangle, need upper left
upperleft = round(xy) - halfwidth;

% need to modify for pixels near edge of images
upper = upperleft(:,2);
left = upperleft(:,1);
lower = upper + default_height;
right = left + default_width;
width = default_width * ones(size(upper));
height = default_height * ones(size(upper));

% check edges for coordinates outside image
[upper,height] = adjust_lo_edge(upper,1,height);
[~,height] = adjust_hi_edge(lower,size(img,1),height);
[left,width] = adjust_lo_edge(left,1,width);
[~,width] = adjust_hi_edge(right,size(img,2),width);

% set width and height to zero when less than default size
iw = find(width<default_width);
ih = find(height<default_height);
idx = unique([iw; ih]);
width(idx) = 0;
height(idx) = 0;

rect = [left upper width height];

%-------------------------------
%
function [coordinates, breadth] = adjust_lo_edge(coordinates,edge,breadth)

indx = find( coordinates<edge );
if ~isempty(indx)
    breadth(indx) = breadth(indx) - abs(coordinates(indx)-edge);
    coordinates(indx) = edge;
end

%-------------------------------
%
function [coordinates, breadth] = adjust_hi_edge(coordinates,edge,breadth)

indx = find( coordinates>edge );
if ~isempty(indx)
    breadth(indx) = breadth(indx) - abs(coordinates(indx)-edge);
    coordinates(indx) = edge;
end

%-------------------------------
%
function [xyinput_in,xybase_in,input,base] = ParseInputs(varargin)

iptchecknargin(4,4,nargin,mfilename);

xyinput_in = varargin{1};
xybase_in = varargin{2};
if size(xyinput_in,2) ~= 2 || size(xybase_in,2) ~= 2
    error(message('images:cpcorr:cpMatrixMustBeMby2'))
end

if size(xyinput_in,1) ~= size(xybase_in,1)
    error(message('images:cpcorr:needSameNumOfControlPoints'))
end

input = varargin{3};
base = varargin{4};
if ndims(input) ~= 2 || ndims(base) ~= 2
    error(message('images:cpcorr:intensityImagesReq'))
end

input = double(input);
base = double(base);

if any(xyinput_in(:)<0.5) || any(xyinput_in(:,1)>size(input,2)+0.5) || ...
   any(xyinput_in(:,2)>size(input,1)+0.5) || ...
   any(xybase_in(:)<0.5) || any(xybase_in(:,1)>size(base,2)+0.5) || ...
   any(xybase_in(:,2)>size(base,1)+0.5)
    error(message('images:cpcorr:cpPointsMustBeInPixCoord'))
end
