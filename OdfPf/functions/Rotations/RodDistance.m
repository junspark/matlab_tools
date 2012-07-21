function dist = RodDistance(pt, ptlist, sym)
% RodDistance - Find angular distance between rotations.
%   
%   USAGE:
%
%   dist = RodDistance(pt, ptlist, sym)
%
%   INPUT:
%
%   pt     is 3 x 1, 
%          a point given in Rodrigues parameters
%   ptlist is 3 x n, 
%          a list of points, also Rodrigues 
%   sym    is 4 x m, 
%          the symmetry group in quaternions
%
%   OUTPUT:
%
%   dist   is 1 x n, 
%          the distance between `pt' and each point in `ptlist'
%
q1   = QuatOfRod(pt);
q2   = QuatOfRod(ptlist);
dist = Misorientation(q1, q2, sym);
%
