function PlotCubicFR(mesh, odf, varargin)
% PLOTCUBICFR - Plot fundamental region for cubic symmetries.
%
%   PlotCubicFR(mesh, odf)
%   PlotCubicFR(mesh, odf, 'param', 'value', ...)
%
%   mesh is a MeshStructure on the fundamental region
%   odf is 1 x n, (where n is the number of nonequivalent
%                 nodes) is the function to plot
%
%   These arguments can be followed by a list of
%   parameter/value pairs which control certain plotting 
%   features.  Options are:
%
%   'ShowMesh'         on|{off}
%                      to show mesh lines on the surface plot
%   'Colormap'         ncolors x 3  (RGB)
%                      to specify the figure's colormap; 
%                      the default colormap is a brightened
%                      version of `jet'
%   'RescaleFigure'    1 x 1 (scalar)
%                      to rescale the figure proportionally
%
%   Notes:
%
%   *)  This function makes a figure containint side-by-side 
%       plots of the cubic fundamental region surface and of 
%       slices through it, with a common colorbar in the middle.
%
%
CheckMesh(mesh, odf)  % check sizes
%
crd = mesh.crd;  % Unpack input structures.
con = mesh.con;
eqv = mesh.eqv;
%
showmesh = 0;         % defaults
%
cmap = [];            % default  colormap set later
bright = 0.4;       
%
plotsurfopts{1} = 'ShowMesh'; % Default for showmesh
plotsurfopts{2} = 'off';
%
figscale  = 1.0;      % figure scale
%
%-------------------- Process options
%
if (nargin >= 3) 
  %
  lenargs = length(varargin);
  nopt    = lenargs/2;
  options = reshape(varargin, 2, nopt);
  %
  for i=1:nopt
    %
    opti = options{1, i};
    vali = options{2, i};
    %
    if (strcmp(opti, 'Colormap'))
      cmap = vali;
      continue;
    end
    %
    if (strcmp(opti, 'ShowMesh')) % passed to PlotSurface
      plotsurfopts{2} = vali;
      continue;
    end
    if (strcmp(opti, 'RescaleFigure'))
      figscale = vali;
      continue;
    end
    %
  end
end
%
%-------------------- Build figure.
%
f = figure;
%
%  Reshape figure to reflect the 1x2 array of subplots.
%
vertscale = 0.6;
%
pos       = get(f, 'Position');
pos(4)    = pos(4)*vertscale;
%
pos = figscale*pos;
%
set(f, 'Position', pos);
%
%-------------------- *** Colormap
%
if (isempty(cmap))
  cmap   = bright + (1 - bright)*jet;
end
%
colormap(cmap);
%
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
clim    = get(a11, 'CLIM');      % save color range for next plot
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
PlotCubicFRPerimeter;
%
%-------------------- *** Subplot 2 (slices)
%
subplot(1, 2, 2)
%
%  The second plot shows the slices and uses the same data
%  as the first plot, so that the colorbar is identical.
%
PlotCubicFRSlices(mesh, odf)
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
