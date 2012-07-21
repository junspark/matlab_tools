function dist = PSphDistance(pt, ptlist)
% PSphDistance - Distance on projective sphere.
%   
%   USAGE:
%
%   dist = PSphDistance(pt, ptlist)
%
%   INPUT:
%
%   pt     is 3 x 1, 
%          a point on the unit sphere (S^2)
%   ptlist is 3 x n, 
%          a list of points on the unit sphere
%
%   OUTPUT:
%
%   dist is 1 x n, 
%        the distance from `pt' to each point in the list
%
%   NOTES:
%
%   *  The distance between two points on the sphere is the angle
%      in radians between the two vectors.  On the projective
%      sphere, antipodal points are considered equal, so the 
%      distance is the minimum of the distances obtained by
%      including the negative of pt as well.
%
%
n = size(ptlist, 2);  % number of points
%
dist2 = [...
    acos(dot(repmat( pt, 1, n), ptlist, 1));...
    acos(dot(repmat(-pt, 1, n), ptlist, 1)) ...
    ];
dist = min(dist2);
