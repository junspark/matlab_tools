function surfs = RectMeshSurfaces(m, varargin)
% RectMeshSurfaces - Surfaces of a rectangular mesh.
%
%   USAGE:
%
%   surfs = RectMeshSurfaces(m, <options>)
%
%   INPUT:
%
%  m is a MeshStructure 
%    It delineates a rectangular region.
%
%   <options> is a sequence of parameter, value pairs
%             Available options are listed below with default values
%             shown in brackets.
%
%    'minx'      {0} | 1
%    'miny'      {0} | 1
%    'minz'      {0} | 1
%    'maxx'      {0} | 1
%    'maxy'      {0} | 1
%    'maxz'      {0} | 1
%
%   OUTPUT:
%
%   surfs is a structure
%         It has components 'minx', etc., for each surface turned
%         on by the options.
%
%   NOTES:
%
%   * Each componenet of surfs is currently just a list of node numbers
%   
MyName = mfilename;
%
%-------------------- Defaults and Keyword Arguments
%
optcell = {...
    'minx',  0, ...
    'miny',  0, ...
    'minz',  0, ...
    'maxx',  0, ...
    'maxy',  0, ...
    'maxz',  0
       };
%
opts = OptArgs(optcell, varargin);
%
%-------------------- *
%
%  Need to find surface elements on boundary and number of 
%  solid element containing it.
%
%
%TBR%switch m.etype.name  % this needs to incorporated with ElementTypeStruct
%TBR% case 'bricks:8'
%TBR%  ELSURFS = [1 4 3 2; 5 6 7 8; 1 5 8 4; 2 3 7 6; 1 2 6 5; 3 4 8 7]';
%TBR% case 'tets:10'
%TBR%  ELSURFS = [1 6 5 4 3 2; 3 4 5 9 10 8; 1 7 10 9 5 6; 1 2 3 8 10 7]';
%TBR% otherwise
%TBR%  error('unexpected element type')
%TBR%end
ELSURFS = m.etype.surfs;
%
sz = size(ELSURFS); nvels = size(m.con, 2); nsels = nvels*sz(2);
%
SURFELS = reshape(m.con(ELSURFS, :), [sz(1), nsels]);
VOLELS  = reshape(repmat(1:nvels, [sz(2), 1]), [1, sz(2)*nvels]);
%
if opts.minx
  tmp = min(m.crd(1, :));
  surfs.minx = find(m.crd(1, :) == tmp);
end
%
if opts.miny
  tmp = min(m.crd(2, :));
  surfs.miny = find(m.crd(2, :) == tmp);
end
%
if opts.minz
  tmp = min(m.crd(3, :));
  surfs.minz = find(m.crd(3, :) == tmp);
end
%
if opts.maxx
  tmp = max(m.crd(1, :));
  surfs.maxx = find(m.crd(1, :) == tmp);
end
%
if opts.maxy
  tmp = max(m.crd(2, :));
  surfs.maxy = find(m.crd(2, :) == tmp);
end
%
if opts.maxz
  tmp = max(m.crd(3, :));
  surfs.maxz = find(m.crd(3, :) == tmp);
end
%
