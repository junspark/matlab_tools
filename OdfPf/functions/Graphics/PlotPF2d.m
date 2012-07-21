function PlotPF2d(smesh, pfdata, varargin)
% PlotPF2d - Plot pole figures in 2D.
%
%   USAGE:
%
%   PlotPF2d(smesh, pfdata, 'param', 'value', ...)
%
%   INPUT:
%
%   smesh   is a MeshStructure,
%           on the sphere or hemisphere
%   pfdata  is a 1 x 2*n cell array,
%           it consists of n pairs, each pair being
%           a set of pole figure nodal point values
%           followed by the pole figure axes (3 x 3
%           orthogonal matrix)
%
%   These arguments can be followed by a list of
%   parameter/value pairs which specify keyword options.
%   Available options include:
%
%   'NumberOfColors'   positive integer
%
%   'RefineLevel'      positive integer
%
%   'TriMeshOpts'      cell array
%                      parameter/value pairs to pass to
%                      trimesh
%   'ColorMap'         ncolors x 3
%                      to specify the figure's colormap
%
%   OUTPUT:  none
%
f = figure;
hold on;
numPfs = length(pfdata)/2;
p = get(f, 'Position');
p(3) = numPfs*p(4);
set(f, 'Position', p);
%
%  Process options.  
%
%  Note, values that are cell arrays need to be placed inside
%  a containing cell array for the call to `struct'.
%
ncolors_dflt = 20;
optkeys = {...
    'NumberOfColors', ncolors_dflt, ...
    'TriMeshOpts',            {{}}, ...  
    'RefineLevel',               3, ...  
    'ColorMap',                 []  ...
	  };
opts = OptArgs(optkeys, varargin);
%
NCOLORS = opts.NumberOfColors;
if (isempty(opts.ColorMap))
  CMAP = jet(NCOLORS);
else
  CMAP = opts.ColorMap;
end
%TRIMESH_OPTS = {'FaceColor', 'interp', 'EdgeColor', 'none', ...
%	       opts.TriMeshOpts{:}};
TRIMESH_OPTS = {'FaceColor', 'flat', 'EdgeColor', 'none', ...
	       opts.TriMeshOpts{:}};
%
xtransl = 4.0;  % space between centers of figures
%
refine_level = opts.RefineLevel;
refined_mesh = RefineMesh(smesh, refine_level);
[refcrd, refcon] = SubdivideSimplex(2, refine_level); % use 3 subdivisions
numel = size(smesh.con, 2);
pts = repmat(refcrd, [1 numel]);
els = repmat(1:numel, [size(refcrd, 2), 1]); els = els(:)';
pts = pts(:, refined_mesh.ord);
els = els(refined_mesh.ord);
%
crd = refined_mesh.crd;
con = refined_mesh.con;
%
ncrd = size(crd, 2);
%
options = {};
%
ldata = length(pfdata);
npf   = ldata/2;
celldata = reshape(pfdata, [2 npf]);
%
allpts = [];
allcon = [];
allpfs = [];
%
for i=1:npf
  %
  pfi = celldata{1, i}; pfi = pfi(:)';
  pfi = EvalMeshFunc(smesh, pfi, els, pts);
  bsi = celldata{2, i};
  %
  
  ptsi = EqualArea(crd, bsi);
  %
  bcrd = bsi'* crd;
  zels = reshape(bcrd(3, con), [3 size(con, 2)]);
  select = (min(zels, [], 1) > -1.0e-6);
  ncon = con(:, select);
  %
  ptsi(1, :) = xtransl*(i-1) + ptsi(1, :);
  %
  allpts = [allpts ptsi];
  allcon = [allcon ncrd*(i-1) + ncon];
  allpfs = [allpfs pfi];
  %
end
%
z = zeros(1, npf*ncrd);
h1 = trimesh(allcon', allpts(1,:), allpts(2, :), z, allpfs, ...
	     TRIMESH_OPTS{:});

hold on;
%
view(0, 90);
colormap(CMAP);
box off;
axis equal;
axis off;
colorbar('EastOutside');
