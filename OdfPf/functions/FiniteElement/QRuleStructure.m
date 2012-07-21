function qrule = QRuleStructure(pts, wts)
% QRuleStructure - Quadrature rule structure.
%   
%   USAGE:
%
%   qrule = QRuleStructure
%   qrule = QRuleStructure(pts, wts)
%
%   INPUT:
%
%   pts is m x n, 
%       a list of n m-dimensional points 
%   wts is an n-vector, 
%       the associated weights
%
%   OUTPUT:
%
%   qrule is a QRuleStructure, 
%         it consists of two fields, points (.pts) and
%         weights (.wts)
%
%   NOTES:
%
%   *  With no arguments, this returns an empty structure.
%
if (nargin == 0)
  pts = [];
  wts = [];
end
%
%  Check for validity.
%
wts = wts(:)';
np = size(pts, 2);
nw = size(wts, 2);
if (np ~= nw)
  error('mismatch in number of points and weights');
end
%
qrule = struct('pts', pts, 'wts', wts);
%
