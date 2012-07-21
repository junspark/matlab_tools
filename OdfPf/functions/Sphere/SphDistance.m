function dist = SphDistance(pt, ptlist)
% SphDistance - Find distance on sphere between a fixed point and others.
%   
%   USAGE:
%
%   dist = SphDistance(pt, ptlist)
%
%   INPUT:
%
%   pt     is 3 x 1, 
%          a point on the unit sphere (S^2)
%   ptlist is 3 x n, 
%          a list of points on the unit sphere (S^2)
%
%   OUTPUT:
%
%   dist is 1 x n, 
%        the distance from `pt' to each point on the list
%
%   NOTES:
%
%   *  The distance is just the angle between the two vectors.
%
n = size(ptlist, 2);  % number of points
%
dist = acos(dot(repmat(pt, 1, n), ptlist, 1));
