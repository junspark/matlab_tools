function PlotDiscretePF(pts, values, basis, range)
% PLOTDISCRETEPF - hey
%   PlotDiscretePF(pts, values, basis, range)
numPts = size(pts, 2);

pts = [EqualChordLengthProj(pts, basis);zeros(1, numPts)];

if nargin > 3
    PlotPoints(pts, values, range)
else
    PlotPoints(pts, values)
end

hold on;
[circ(:, 1), circ(:, 2)] = pol2cart([0:2*pi/999:2*pi]', sqrt(2));
plot(circ(:, 1), circ(:, 2), 'k-', 'LineWidth', 1.5)
%
view(0, 90);
%
axis equal
%
hold off