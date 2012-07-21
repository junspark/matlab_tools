function xyz = XYZOfThetaPhi(thetaphi)
% XYZOfThetaPhi - Map spherical coordinates to sphere.
%   
%   USAGE:
%
%   xyz = XYZOfThetaPhi(thetaphi)
%
%   INPUT:
%
%   thetaphi is 2 x n, 
%            the spherical coordinates for a list of n points; 
%            theta is the angle that the projection onto x-y plane 
%            makes with the x-axis, and phi is the angle with z-axis
%
%   OUTPUT:
%
%   xyz is 3 x n, 
%       the Cartesian coordinates of the points described by (theta, phi)
%
%   NOTES: 
%
%   *  The matlab builtin `sph2cart' could also be used, but 
%      the convention for phi is different (90 degrees minus thisone).
%
theta = thetaphi(1, :);
phi   = thetaphi(2, :);
%
ct = cos(theta); 
st = sin(theta);
%
sp = sin(phi);
%
xyz = [sp .* ct ; sp .* st; cos(phi)];
