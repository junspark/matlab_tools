function gij = MetricGij(diff)
% MetricGij - Compute components of metric from differential.
%   
%   USAGE:
%
%   gij = MetricGij(diff)
%
%   INPUT:
%
%   diff is m x n x l,
%        the array of n tangent vectors of dimension m at each 
%        of l points
%
%   OUTPUT:
%
%   gij is n x n x l, 
%       the metric components (dot(ti, tj)) at each of the l points
%
gij = MultMatArray(permute(diff, [2 1 3]), diff);
