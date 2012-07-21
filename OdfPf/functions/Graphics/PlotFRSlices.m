function PlotFRSlices(mesh, odf, type, planes)
% PlotFRSlices - Plot slices of the fundamental region.
%   
%   VERSION:  $Id: PlotFRSlices.m 10 2008-04-10 15:27:09Z boyce $
%
%   USAGE:
%
%   PlotFRSlices(mesh, odf, type)
%   PlotFRSlices(mesh, odf, type, planes)
%
%   INPUT:
%
%   mesh is a MeshStructure,
%        on the given fundamental region
%   odf  is a vector,
%        the odf values at the independent nodal points of `mesh'
%   type is a string, (optional, default: 'cubic')
%        indicating the symmetry type; possible values are
%        'cubic' | 'hexagonal'; the slices to display
%        are based on the type
%   planes is an array of structures, (optional)
%        each structure is expected to have components 'Point'
%        and 'Normal' characterizing a plane;  if specified,
%        this overrides the `type' input
%
%   OUTPUT:  none
%
if (nargin < 3)
  type = 'cubic';
end
%
if (nargin < 4)
  planes = SymPlanes(type);
end
%
if isempty(planes)
  planes = SymPlanes('cubic');
end
%

edgeopts = {'EdgeColor', 'black'};

for p=planes
  [sm, e, ecrd] = SliceMesh(mesh, p.Point, p.Normal);

  if length(odf) < size(mesh.crd, 2)
    odf = ToAllNodes(odf, mesh.eqv);
  end
  eodf      = odf(mesh.con(:, e));
  slice_odf = dot(ecrd, eodf, 1);
  %
  PlotSurface(sm, slice_odf);
  %
end
%
%--------------------Private Functions----------------------------------
%
function planes = SymPlanes(type)
% CUBPLANES - 
%   
Z  = [0; 0; 0];
%
if (strcmpi(type, 'cubic'))
  %
  E  = eye(3);
  E1 = E(:, 1);
  E2 = E(:, 2);
  E3 = E(:, 3);
  planes = struct('Point', Z, 'Normal', {E1, E2, E3});
  %
elseif (strcmpi(type, 'hexagonal'))
  %
  n_z  = [0 0 1]';
  theta = (0:5)*2*pi/5;
  n_xy = [cos(theta); sin(theta); zeros(1, 6)];
  planes = struct('Point', Z, 'Normal', num2cell([n_z, n_xy], 1));
  %
else
  %
  %  Default planes to display.
  %
  E  = eye(3);
  E1 = E(:, 1);
  E2 = E(:, 2);
  E3 = E(:, 3);
  planes = struct('Point', Z, 'Normal', {E1, E2, E3});
  %
end
