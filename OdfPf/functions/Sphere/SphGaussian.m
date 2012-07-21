function fsm = SphGaussian(cen, pts, stdev)
% SphGaussian - Gaussian distribution on angular distance.
%   
%   USAGE:
%
%   fsm = SphGaussian(cen, pts, stdev)
%
%   INPUT:
%
%   cen   is 3 x 1, 
%         the center of the distribution
%   pts   is 3 x n, 
%         a list of points on the sphere at which to
%         evaluate the result
%   stdev is a scalar, 
%         the standard deviation
%
%   OUTPUT:
%
%   fsm is 1 x n, 
%       the list of values at each point of `pts'
%
%   NOTES:
%
%   *  This returns the values of a 1D Gaussian applied
%      to spherical distance from the center point
%   *  The result is not normalized to have unit integral
%      over the sphere.
%
twosigsq = 2*(stdev^2);
theta    = SphDistance(cen, pts);
minusthetasq = -theta.*theta;
%
fsm = (1/(stdev*sqrt(2*pi)))*exp(minusthetasq/twosigsq);
