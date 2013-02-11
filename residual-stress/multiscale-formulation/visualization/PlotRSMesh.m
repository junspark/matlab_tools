function f = PlotRSMesh(f, np, x, y, z, c, varargin)
% PlotRSMesh - plot continuum length scale data from multiscale residual
% stress formulation on continuum mesh
%
%   USAGE:
%
%   PlotRSMesh(f, np, x, y, z, c)
%   PlotRSMesh(f, np, x, y, z, c, varargin)
%
%   INPUT:
%
%   f
%       figure object
%
%   np
%       mesh connectivity (n x 8)
%       supports 8 noded brick only at the moment
%
%   x
%       mesh x-coordinates (n x 1)
%
%   y
%       mesh y-coordinates (n x 1)
%
%   z
%       mesh z-coordinates (n x 1)
%
%   c
%       field data over mesh (n x 1)
%
%   These arguments can be followed by a list of
%   parameter/value pairs which control certain plotting
%   features.  Options are:
%
%   'Title'            title of the figure (default is 'MESH')
%   'ViewAngle'        angle to view the pole figure from
%                      default is [45 30]
%   'DataRange'        range of the Data to determine the colors
%                      default is [-600 1200]
%   'ShowMesh'         on|{off}
%                      to display the mesh
%   'Colorbar'         on|{off}
%                      to display the colorbar
%
%   OUTPUT:
%
%   f
%       figure object
%

% default options
opts    = {...
    'Title', 'MESH', ...
    'ViewAngle', [45 30], ...
    'DataRange', [-600 1200], ...
    'ShowMesh', 'off', ...
    'Colorbar', 'off', ...
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
    zz  = [z(n1) z(n2) z(n3) z(n4)...
        z(n5) z(n6) z(n7) z(n8)];
    cc  = [c(n1) c(n2) c(n3) c(n4)...
        c(n5) c(n6) c(n7) c(n8)];
    
    if strcmp(opts.ShowMesh, 'on')
        line([xx(1:4) xx(1)],[yy(1:4) yy(1)],[zz(1:4) zz(1)],'Color','k');
        line([xx(5:8) xx(5)],[yy(5:8) yy(5)],[zz(5:8) zz(5)],'Color','k');
        line([xx(1) xx(5)],[yy(1) yy(5)],[zz(1) zz(5)],'Color','k');
        line([xx(2) xx(6)],[yy(2) yy(6)],[zz(2) zz(6)],'Color','k');
        line([xx(3) xx(7)],[yy(3) yy(7)],[zz(3) zz(7)],'Color','k');
        line([xx(4) xx(8)],[yy(4) yy(8)],[zz(4) zz(8)],'Color','k');
        plot3(xx, yy, zz, 'k.')
        hold on
    end
    
    % FACE 1
    xp  = xx(1:4); yp  = yy(1:4); zp  = zz(1:4); cp  = cc(1:4);
    p   = patch(xp, yp, zp, cp, ...
        'EdgeColor', 'none', 'LineStyle', '-');
    hold on
    
    % FACE 2
    xp  = xx(5:8); yp  = yy(5:8); zp  = zz(5:8); cp  = cc(5:8);
    p   = patch(xp, yp, zp, cp, ...
        'EdgeColor', 'none', 'LineStyle', '-');
    
    % FACE 3
    xp  = [xx(1) xx(5) xx(8) xx(4)];
    yp  = [yy(1) yy(5) yy(8) yy(4)];
    zp  = [zz(1) zz(5) zz(8) zz(4)];
    cp  = [cc(1) cc(5) cc(8) cc(4)];
    p   = patch(xp, yp, zp, cp, ...
        'EdgeColor', 'none', 'LineStyle', '-');
    
    % FACE 4
    xp  = [xx(2) xx(3) xx(7) xx(6)];
    yp  = [yy(2) yy(3) yy(7) yy(6)];
    zp  = [zz(2) zz(3) zz(7) zz(6)];
    cp  = [cc(2) cc(3) cc(7) cc(6)];
    p   = patch(xp, yp, zp, cp, ...
        'EdgeColor', 'none', 'LineStyle', '-');
    
    % FACE 5
    xp  = [xx(2) xx(6) xx(5) xx(1)];
    yp  = [yy(2) yy(6) yy(5) yy(1)];
    zp  = [zz(2) zz(6) zz(5) zz(1)];
    cp  = [cc(2) cc(6) cc(5) cc(1)];
    p   = patch(xp, yp, zp, cp, ...
        'EdgeColor', 'none', 'LineStyle', '-');
    
    % FACE 6
    xp  = [xx(3) xx(7) xx(8) xx(4)];
    yp  = [yy(3) yy(7) yy(8) yy(4)];
    zp  = [zz(3) zz(7) zz(8) zz(4)];
    cp  = [cc(3) cc(7) cc(8) cc(4)];
    p   = patch(xp, yp, zp, cp, ...
        'EdgeColor', 'none', 'LineStyle', '-');
end
axis equal tight off

title(opts.Title, ...
    'FontSize', 12, 'FontWeight', 'bold')

view(opts.ViewAngle);

caxis(opts.DataRange);

if strcmp(opts.Colorbar, 'on')
    colorbar('location', 'EastOutside')
end
