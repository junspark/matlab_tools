function PlotSphere(smesh, pf, varargin)
% PlotSphere - Plot a function (pole figure) on the Sphere
%
%   USAGE:
%
%   PlotSphere(smesh, pf)
%   PlotSphere(smesh, pf, 'param', 'value', ...)
%
%   INPUT:
%
%   smesh is a MeshStructure,
%         on the sphere (the usual one)
%   pf    is a vector, 
%         the nodal point values of thefunction (pole figure) to plot
%
%   These arguments can be followed by a list of
%   parameter/value pairs which control certain plotting 
%   features.  Options are:
%
%   'ShowMesh'         on|{off}
%                      to show mesh lines
%   'ShowQuadrants'    {on}|off
%                      to show quadrant divisions on
%                      the sphere
%   'UpperHemisphere'  on|{off}
%                      to limit quadrant divisions to 
%                      only the upper hemisphere; this only
%                      applies if 'ShowQuadrants' is in
%                      effect
%   'Colormap'         ncolors x 3  (RGB)
%                      to specify the figure's colormap; 
%                      the default colormap is a brightened
%                      version of `jet'
%   'ShowColorBar'     on|{off}
%                      to display a colorbar or not
%
%   OUTPUT:  none
%
%   NOTES:
%
%   *  `hold on' is in effect following this routine
%
CheckMesh(smesh, pf)  % check array sizes
%
%
%--------------------Defaults and Options-------------------------------
%
optcell = {...
    'ShowMesh',         'off',    ...
    'ShowQuadrants',     'on',    ...
    'UpperHemisphere',  'off',    ...
    'ShowColorBar',      'on',    ...
    'Colormap',            [],    ...
    'BrightenColormap',     0,    ...
    'NumberOfColors',      64     ...
    };
%
opts = OptArgs(optcell, varargin);
%

%
showmesh = OnOrOff(opts.ShowMesh);
showquad = OnOrOff(opts.ShowQuadrants);
upperhem = OnOrOff(opts.UpperHemisphere);
showcbar = OnOrOff(opts.ShowColorBar);
%
%-------------------- *** Colormap
%
f = figure;
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
colormap(cmap);
%
%  Surface options.
%
plotsurfopts = {'ShowMesh', opts.ShowMesh};
%
%-------------------- Plot surface
%
PlotSurface(smesh, pf, plotsurfopts{:})   
%
hold on;
%
%-------------------- Show quadrant divisions.
%
if (showquad)
  PlotSphereQuadrants(100, upperhem);
end
%
%-------------------- Colorbar
%
if (showcbar == 1)
  colorbar;
end

