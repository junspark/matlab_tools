function PlotSurface(smesh, sfun, varargin)
% PlotSurface - Plot a function on surface in 3d.
%
%   USAGE:
%
%   PlotSurface(smesh, sfun)
%   PlotSurface(smesh, sfun, 'param', 'value', ...)
%
%   INPUT:
%
%   smesh is a MeshStructure,
%         on a 2D surface in 3D
%   sfun  is a vector, 
%         the values on the nodes of the surface mesh
%
%   These arguments can be followed by a list of
%   parameter/value pairs which control certain plotting 
%   features.  Options are:
%   
%   'ShowMesh'    on|{off}
%                 to show the edges of the mesh
%
%   OUTPUT:  none
%
%   NOTES:
%
%   *  This routine calls `trimesh' with interpolated face color.
%
optcell = {...
    'ShowMesh',  'off' ...
    };
%
opts = OptArgs(optcell, varargin);
%
showmesh = OnOrOff(opts.ShowMesh);
% 
scrd = smesh.crd;
scon = smesh.con;
seqv = smesh.eqv;
%
CheckMesh(smesh, sfun)    %  Check consistency of function data.
%
sfunall = ToAllNodes(sfun, seqv);
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
    if (strcmp(opti, 'ShowMesh'))
      showmesh = OnOrOff(vali);
      continue;
    end
  end
end
%
%-------------- Plot function on mesh.
%
if (showmesh == 1)
  edgeopts = {'EdgeColor', 'black'};
else
  edgeopts = {'EdgeColor', 'none', 'LineStyle', '-'};
end
%
h1 = trimesh(scon', scrd(1,:), scrd(2,:), scrd(3,:), sfunall, ...
	      'FaceColor', 'interp', ...
	      edgeopts{:});
hold on;
%
%-------------- Show mesh lines if requested.
%
%if (showmesh == 1)
%  %
%  drawscale = 1.001;
%  %
%  h2 = patch('Faces',     scon',                  ...
%	      'Vertices',  drawscale*scrd',        ...
%	      'FaceColor', 'none',                 ...
%	      'EdgeColor', 'black',                ...
%	      'LineWidth',     2.0                 ...
%	      );
%end
%
%-------------- Axes properties.
%
axis equal
axis off
%
%set(gca, 'Color', 'none', 'Visible', 'off');
% set(gca, 'CLim', [0 5]);
