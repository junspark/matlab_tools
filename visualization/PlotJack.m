function figobj = PlotJack(input_matrix, varargin)
% PlotJack - plot a jack for a symmetric 3 x 3 matrix
%
%   USAGE:
%
%   PlotJack(input_matrix)
%   PlotJack(input_matrix, varargin)
%
%   INPUT:
%
%   input_matrix
%       symmetric 3 x 3 matrix
%
%   This argument can be followed by a list of
%   parameter/value pairs which control certain plotting
%   features.  Options are:
%
%   'ViewAngle'        angle to view the jack
%                      default is [45 30]
%   'DataRange'        range of the Data to determine the colors
%                      default is [minimum(eigen_value(input_matrix)) maximum(eigen_value(input_matrix))]
%   'JackLocation'     location of the jack's origin
%                      default is [0 0 0]
%   'xaxis'            label for the x-axis
%                      default is 'x'
%   'yaxis'            label for the y-axis
%                      default is 'y'
%   'zaxis'            label for the z-axis
%                      default is 'z'
%
%   OUTPUT:
%
%   figobj
%       figure object of the jack plot
%

issym = all(all(input_matrix == input_matrix.'));
if ~issym
    error('input matrix is not symmetric')
end

[eigen_direction, eigen_value]  = eig(input_matrix);

eigen_value                 = diag(eigen_value);
[eigen_value, col_index]    = sort(eigen_value,1,'ascend');
eigen_direction             = eigen_direction(:,col_index);

DataRange   = [eigen_value(1) eigen_value(end)];

% default options
optcell = {...
    'ViewAngle', [45 30], ...
    'DataRange', DataRange, ...
    'JackLocation', [0 0 0]', ...
    'xaxis', 'x', ...
    'yaxis', 'y', ...
    'zaxis', 'z', ...
    };

% update option
opts    = OptArgs(optcell, varargin);

ncmap   = 256;
cmap    = jet(ncmap);
dData   = opts.DataRange(2) - opts.DataRange(1);

Ox  = opts.JackLocation(1);
Oy  = opts.JackLocation(2);
Oz  = opts.JackLocation(3);

figobj  = figure;
hold on
for i = 1:1:3
    QuiverColor = round((ncmap-1)*(eigen_value(i) - DataRange(1))./dData) + 1;
    
    Vxp = 0.5*eigen_direction(1,i);
    Vyp = 0.5*eigen_direction(2,i);
    Vzp = 0.5*eigen_direction(3,i);
    
    Vxn = -0.5*eigen_direction(1,i);
    Vyn = -0.5*eigen_direction(2,i);
    Vzn = -0.5*eigen_direction(3,i);
    
    quiver3(Ox, Oy, Oz, Vxp, Vyp, Vzp, ...
        'Color', cmap(QuiverColor, :), 'LineWidth', 5);
    quiver3(Ox, Oy, Oz, Vxn, Vyn, Vzn, ...
        'Color', cmap(QuiverColor, :), 'LineWidth', 5);
end

colorbar('Location', 'EastOutside', 'FontSize', 16, 'FontWeight', 'bold')
caxis(opts.DataRange);
axis square equal
grid on
quiver3(0, 0, 0, 1, 0, 0, ...
    'Color', 'k', ...
    'LineWidth', 2)
text(1.1, 0, 0, opts.xaxis, ...
    'FontSize', 16, 'FontWeight', 'bold')
quiver3(0, 0, 0, 0, 1, 0, ...
    'Color', 'k', ...
    'LineWidth', 2)
text(0, 1.1, 0, opts.yaxis, ...
    'FontSize', 16, 'FontWeight', 'bold')
quiver3(0, 0, 0, 0, 0, 1, ...
    'Color', 'k', ...
    'LineWidth', 2)
text(0, 0, 1.1, opts.zaxis, ...
    'FontSize', 16, 'FontWeight', 'bold')
view(opts.ViewAngle);