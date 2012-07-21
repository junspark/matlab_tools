function PlotSphereQuadrants(ndiv, upper)
% PlotSphereQuadrants - Outline quadrants on unit sphere
%
%   USAGE:
%
%   PlotSphereQuadrants(ndiv)
%   PlotSphereQuadrants(ndiv, upper)
%   
%   INPUT:
%
%   ndiv  is a positive integer, 
%         the number of points per circle
%   upper is a positive integer, (optional, default=0) 
%         if nonzero, only the upper half is drawn
%
%   OUTPUT:  none
%
if (nargin < 2)    % process arguments
  upper = 0;
end
%
%-------------------- Draw quadrants.
%
scale = 1.001;
%
t = (pi/ndiv)*(0:ndiv);
%
c = scale*cos(t)';
s = scale*sin(t)';
z = zeros(size(c));
%
x1 = [c -c  c  z]; 
x2 = [s -s  z  c];
x3 = [z  z  s  s];
%
hold on;
%
plot3(x1, x2, x3, 'k-', ...
      'LineWidth', 1.0);
%
if (~upper)
  %
  plot3(-x1, -x2, -x3, 'k--', ...
	'LineWidth', 1.0);
  %
end
%
