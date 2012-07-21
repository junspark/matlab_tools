function stats = MisorientationStats(misorient, locations, wts)
% MisorientationStats - Misorientation correlation statistics.
%
%   USAGE:
%
%   stats = MisorientationStats(misorient, locations)
%   stats = MisorientationStats(misorient, locations, wts)
%
%   INPUT:
%
%   misorient is 4 x n, 
%             a list of misorientation quaternions,
%             assumed to have been derived from properly clustered 
%             orientation data
%   locations is d x n, (d <= 3) 
%             a list of spatial locations corresponding to the 
%             misorientations
%   wts       is 1 x n, (optional)
%             a list of weights; if not specified, uniform weights are used
%
%   OUTPUT:
%
%   stats is a structure with five components:
%
%         W     is a 3 x 3 matrix (A in Barton paper)
%         X     is a d x d matrix (M in Barton paper)
%         WX    is a 3 x d matrix (cross-correlation
%                    of normalized variables; X in
%                    Barton paper)
%         wi    is 3 x n, the unnormalized axial vectors
%         xi    is d x n, the unnormalized spatial directions
%                         from the centroid
%
%   REFERENCE:  
%
%   "A Methodology for Determining Average Lattice Orientation and 
%   Its Application to the  Characterization of Grain Substructure",
%
%   Nathan R. Barton and Paul R. Dawson,
%
%   Metallurgical and Materials Transactions A,
%   Volume 32A, August 2001, pp. 1967--1975
%   
[d, n] = size(locations);
%
if nargin < 3
  wts = repmat(1/n, [3 n]);
else
  wts = repmat(wts(:)', [3 1]); % make into row vector
end
%
ang   = 2*acos(misorient(1,:));
limit = (ang < eps);
wsc(~limit) = (ang/sin(ang/2));
wsc( limit) = 2.0;
wi  = misorient(2:4, :).* repmat(wsc, [3 1]);
%
cen = sum(locations, 2)/n;
xi  = locations - repmat(cen, [1 n]);
%
Winv = sum(RankOneMatrix(wi.*wts, wi), 3);
Xinv = sum(RankOneMatrix(xi.*wts, xi), 3);
%
W = inv(Winv); Whalf = sqrtm(W);
X = inv(Xinv); Xhalf = sqrtm(X);
%
wibar = (Whalf*wi);
xibar = (Xhalf*xi);
%
WX = sum(RankOneMatrix(wibar.*wts, xibar), 3);
%
stats.W  = W;
stats.X  = X;
stats.WX = WX;
stats.wi = wi;
stats.xi = xi;
