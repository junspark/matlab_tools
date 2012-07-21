function PlotFR(mesh, odf, varargin)
% PlotFR - Plot a function on a fundamental region.
%
%   USAGE:
%
%   PlotFR(mesh, odf)
%   PlotFR(mesh, odf 'param', 'value', ...)
%
%   INPUT:
%
%   mesh is a MeshStructure:
%           describing the fundamental region to be plotted
%   odf  is a vector of reals:,
%           it contains nodal point values of the function to be
%           plotted;  the values are expected only for the
%           independent nodes of the mesh
%
%   These arguments can be followed by a list of
%   parameter/value pairs which control certain plotting
%   features.  Options are:
%
%   'Symmetries'       string indicating symmetries of the FR
%                      {'cubic'} | 'hexagonal' | 'orthorhombic'
%   'ShowMesh'         on | {off}
%                      to show mesh lines on the surface plot
%   'CameraPosition'   sets the camera position for the plots 
%                      using campos
%   'Colormap'         ncolors x 3  (RGB)
%                      to specify the figure's colormap;
%                      the default colormap is a brightened
%                      version of `jet'
%   'Clim'             clim min and max
%   'BrightenColormap' scalar between 0 and 1 (default=0)
%                      factor for brightening the colormap;
%                      a value of 0 does nothing, a value
%                      of 1 makes everything white
%   'NumberOfColors'   positive integer (default = 64)
%                      number of colors to use in default colormap
%
%   OUTPUT:  none
%
%   NOTES:
%
%   *  This function makes a figure containing side-by-side
%      plots of the surface of the given fundamental region and of
%      slices through it, with a common colorbar in the middle.
%
%

% campos([3.7, 5.5, 2.6]) %%% ANDY CAM POS
CheckMesh(mesh, odf)  % check sizes
%
crd = mesh.crd;  % Unpack input structures.
con = mesh.con;
eqv = mesh.eqv;
%
%odf = ToAllNodes(odf, eqv);
%
%--------------------Defaults and Options-------------------------------
%
optcell = {...
    'Symmetries', 'cubic', ...
    'ShowMesh',  'off',    ...
    'CameraPosition', [-3.7824 -4.9293 3.5872], ...
    'Colormap',  [],       ...
    'Clim',  [],       ...
    'BrightenColormap', 0, ...
    'NumberOfColors',  64  ...
    };
%
opts = OptArgs(optcell, varargin);
%
symtype = opts.Symmetries;
%
%  Colormap.
%
cmap = opts.Colormap;
if (isempty(cmap))  % default colormap
    cmap = jet(opts.NumberOfColors);
end
bright = opts.BrightenColormap;
cmap   = bright + (1 - bright)*cmap;
%
%  Surface options.
%
plotsurfopts = {'ShowMesh', opts.ShowMesh};
%
%-------------------- Build figure.
%
f = figure;
set(f, 'PaperPosition', [0 0 4.444 2]);
%
colormap(cmap);
%
figscale  = 1.0;      % figure scale
%
%  Reshape figure to reflect the 1x2 array of subplots.
%
vertscale = 0.6;
%
pos       = get(f, 'Position');
pos(4)    = pos(4)*vertscale;
%
set(f, 'Position', pos);
%-------------------- *** Subplot 1 (surface)
%
subplot(1, 2, 1)
%
%  The first plot shows the surface, and the range of data
%  is complete since the surface plot uses the same nodal point
%  array, but only the surface elements.
%
[faces, multiplicity] = MeshFaces(con);
scon = faces(:, (multiplicity == 1));
surfmesh = MeshStructure(crd, scon, eqv);
%
PlotSurface(surfmesh, odf, plotsurfopts{:});
%
%  Axes.
%
axis off
%
a11     = gca;
clim    = opts.Clim;
if (isempty(clim))  
    clim    = get(a11, 'CLIM');
else
    set(a11, 'CLIM', clim)
end

%
%  Colorbar.
%
%  This positions the colorbar in the middle of the figure.
%  The vertical size is 90% that of the axes, while the default
%  width remains the same.
%
%  Position = [left bottom width height]
%
cb      = colorbar;
pos_a11 = get(a11, 'Position');  % save position after creation of colorbar
poscb   = get(cb, 'Position');
%
poscb(4) = 0.9*pos_a11(4);
poscb(1) = 0.5 - 0.5*poscb(3);
poscb(2) = 0.5 - 0.5*poscb(4);
%
set(cb, 'Position', poscb);
%
PlotFRPerimeter(symtype);
campos(opts.CameraPosition)
%
%-------------------- *** Subplot 2 (slices)
%
subplot(1, 2, 2)
%
%  The second plot shows the slices and uses the same data
%  as the first plot, so that the colorbar is identical.
%
PlotFRSlices(mesh, ToAllNodes(odf, eqv), symtype)
campos(opts.CameraPosition)
% 
%  Adjust axes.
%
a12 = gca;
%
%  Make same size as other axes.
%
pos_a12    = pos_a11;
pos_a12(1) = 1 - pos_a12(1) - pos_a12(3);
set(a12, 'CLIM', clim, 'Position', pos_a12);
