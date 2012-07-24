function [f] = PlotMesh2D(f, np, x, y, c, varargin)

% default options
opts    = {...
    'Title', 'MESH', ...
    'ViewAngle', [90 0], ...
    'DataRange', [-500 500], ...
    'ShowMesh', 'off', ...
    };

% update option
opts    = OptArgs(opts, varargin);

numel   = size(np,1);
for i = 1:1:numel
    n1  = np(i,1);
    n2  = np(i,2);
    n3  = np(i,3);
    n4  = np(i,4);
    n5  = np(i,5);
    n6  = np(i,6);
    n7  = np(i,7);
    n8  = np(i,8);
    
    xx  = [x(n1) x(n2) x(n3) x(n4)...
        x(n5) x(n6) x(n7) x(n8)];
    yy  = [y(n1) y(n2) y(n3) y(n4)...
        y(n5) y(n6) y(n7) y(n8)];
    cc  = [c(n1) c(n2) c(n3) c(n4)...
        c(n5) c(n6) c(n7) c(n8)];
    
    if strcmp(opts.ShowMesh, 'on')
        line([xx xx(1)],[yy yy(1)],'Color','k');
        plot3(xx, yy, zz, 'k.')
        hold on
    end
    
    % FACE 1
    p   = patch(xx, yy, cc, ...
        'EdgeColor', 'none', 'LineStyle', '-');
    hold on
end
axis equal tight