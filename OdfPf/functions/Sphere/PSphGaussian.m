function fsm = PSphGaussian(center, pts, stdev)
% PSphGaussian - Gaussian distribution for smoothing on projective sphere.
%   
%   USAGE:
%
%   fsm = PSphGaussian(center, pts, stdev)
%
%   INPUT:
%
%   center is 3 x 1, 
%          the center of the distribution
%   pts    is 3 x n, 
%          a list of points on the sphere; antipodal
%          points are considered equal
%   stdev  is 1 x 1, 
%          the (1D) standard deviation
%
%   OUTPUT:
%
%   fsm is 1 x n, 
%       the list of values at each point of pts
%
%   Notes:  
%
%   *  The result is not normalized, so this may have to be
%      done after the fact.
%   *  The distribution is a 1-D normal distribution applied
%      to the distance function on the projective sphere.
%   *  The actual scaling factor to give unit integral over 
%      the projective sphere is not computed; the result
%      is not normalized.
%
twosigsq = 2*(stdev^2);
theta    = PSphDistance(center, pts);
minusthetasq  = -theta.*theta;
%
fsm = (1/(stdev*sqrt(2*pi)))*exp(minusthetasq/twosigsq);
