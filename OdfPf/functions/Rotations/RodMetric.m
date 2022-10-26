function metric = RodMetric(rod)
% RodMetric - find volume integration factor due to metric
%   
%   USAGE:
%
%   metric = RodMetric(rod)
%   
%   INPUT:
%
%   rod is 3 x n, 
%       an array of 3-vectors (Rodrigues parameters)
%
%   OUTPUT:
%
%   metric is 1 x n, 
%          the metric at each point of `rod'
%
metric = (1./(1+dot(rod, rod))).^2;
