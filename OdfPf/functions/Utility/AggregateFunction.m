function aggf = AggregateFunction(pts, agg, wts, PointFun, varargin)
% AggregateFunction - Create a function from an aggregate of points.
%   
%   USAGE:
%
%   aggf = AggregateFunction(pts, agg, wts, @PointFun)
%   aggf = AggregateFunction(pts, agg, wts, @PointFun, pfarg1, ...)
%
%   INPUT:
%
%   pts is d x n, 
%       the set of points on which to evaluate the aggregate function
%   agg is d x m, 
%       a collection of points (the aggregate)
%   wts is 1 x m, 
%       the weighting associated with points in `agg'
%   PointFun is a function handle, 
%       the function which evaluates the distribution associated 
%       with each point of the aggregate;
%       the required interface to PointFun is:
%
%       PointFun(center, points [, args])
%                center is d x 1, the center of the distribution
%                points is d x n, a list of points to evaluate
%                args are function-specific arguments
%
%   Remaining arguments are passed to PointFun.
%
%   OUTPUT:
%
%   aggf is 1 x n, 
%        the values of the aggregate function at each point in `pts';
%
%   NOTES:
%
%   * Each point in the aggregate is the center of a distribution
%     over the whole space, given by PointFun; all of these 
%     distributions are superposed to give the resulting 
%     aggregate function, which is then evaluated at the 
%     specified point.
%
[d, n]      = size(pts);
[dcheck, m] = size(agg);
wtscheck = length(wts);
%
if (m ~= wtscheck) 
  error('dimension mismatch:  wts and agg (length)')
end
%
if (d ~= dcheck) 
  error('dimension mismatch:  pts and agg (first dimension)')
end
%
aggf = zeros(1, n);
%
for i=1:m
  aggf = aggf + wts(i) * feval(PointFun, agg(:, i), pts, varargin{:});
end
