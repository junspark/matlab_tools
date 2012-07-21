function qrule = LoadQuadrature(qname)
% LoadQuadrature - Load quadrature data rule.
%   
%   USAGE:
%
%   qrule = LoadQuadrature(qname)
%
%   INPUT:
%
%   qname is a string, 
%         the basename of the quadrature data files
%
%   OUTPUT:
%
%   qrule is a QRuleStructure, 
%         it consists of the quadrature point locations and weights
%
%   NOTES:
%
%   *  It is expected that the quadrature rules are for simplices,
%      and the last barycentric coordinate is appended to the file
%      read from the data.
%
%   *  Expected suffixes are .qps for the location and .qpw for
%      the weights.
%
try
  pts = load([qname '.qps'])';
  wts = load([qname '.qpw'])';
catch
  error('failed to open quadrature data files');
end
%
n   = size(pts, 2);
pts = [pts; ones(1, n) - sum(pts, 1)]; % last barycentric coordinate
%
qrule = QRuleStructure(pts, wts);
