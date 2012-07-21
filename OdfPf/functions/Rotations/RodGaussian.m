function gauss = RodGaussian(cen, pts, stdev, sym)
% RODGAUSSIAN - Gaussian distribution on angular distance.
%   
%   USAGE:
%
%   gauss = RodGaussian(cen, pts, stdev, sym)
%
%   INPUT:
%
%   cen   is 3 x 1, 
%         the center of the distribution (in Rodrigues parameters)
%   pts   is 3 x n, 
%         a list of points (Rodrigues parameters)
%   stdev is 1 x 1, 
%         the standard deviation of the distribution
%   sym   is 4 x k, 
%         the symmetry group (quaternions)
%
%   OUTPUT:
%
%   gauss is 1 x n, 
%         the list of values at each input point
%
%   NOTES:
%
%   *  This returns the values of a (not normalized) 1D Gaussian 
%      applied to angular distance from the center point 
%   *  The result is not normalized to have unit integral.
%
twosigsq     = 2*(stdev^2);
theta        = RodDistance(cen, pts, sym);
minusthetasq = -theta.*theta;
%
%gauss = (1/(stdev*sqrt(2*pi)))*exp(minusthetasq/twosigsq);
gauss = exp(minusthetasq/twosigsq);
