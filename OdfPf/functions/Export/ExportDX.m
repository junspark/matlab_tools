function header = ExportDX(fname, mesh, Ndata, Edata, varargin)
% ExportDX - Write a DX field in DX native format.
%   
%   USAGE:
%
%   header = ExportDX(fname, mesh)
%   header = ExportDX(fname, mesh, Ndata)
%   header = ExportDX(fname, mesh, Ndata, Edata)
%
%   INPUT:
%
%   fname is a string, 
%         the basename of the DX header file, as well as the 
%         basename of the various ASCII output files, and also
%         the name of the resulting DX field
%
%   mesh  is a MeshStructure
%
%   Ndata is a cell array,
%         it contains position-dependent data; it has the form
%         {name1, array1, name2, array2, ...}, where the names
%         are strings with no spaces and the arrays are real-valued
%         data arrays;  Ndata can be empty.
%         
%   Edata is a cell array,
%         it contains connection-dependent data and has the same
%         form as Ndata
%
%   These arguments can be followed by a list of
%   parameter/value pairs which specify keyword options.
%   Available options include:
%
%   'ElementType'     string specifying element type; 
%                     A default value is determined from the
%                     connectivity array, but for four nodes 
%                     per element, 'tetrahedra' is chose over 
%                     'quads'; possible values are:
%                     'lines'|'triangles'|'quads'|'tetrahedra'|'cubes'
%
%   OUTPUT:
%
%   header is a cell array of strings,
%          it contains the DX header file content
%
%   NOTES:
%
%   *  ExportDX writes the DX field to fname.dx, and the field 
%      object is named fname.  Each array in Ndata or Edata is
%      written to a separate file (fname-arayname.dat) and is 
%      added as a component of the DX field.
%
%   *  Currently, this routine can only handle element types
%      of 'lines', 'triangles', 'tetrahedra' or 'cubes'.
%
if (nargin < 2)      % check arguments
  error('Need at least two arguments: fname, mesh')
elseif (nargin == 2)
  Ndata = {};
  Edata = {};
elseif (nargin == 3)
  Edata = {};
end
%
%-------------------- Defaults and Keyword Arguments
%
optcell = {...
    'ElementType',  '' ...
	  };
%
opts = OptArgs(optcell, varargin);
%
%-------------------- Defaults and Keyword Arguments
%
%
crd = mesh.crd';
con = mesh.con'; con = con - 1; %  DX is 0-based
eqv = mesh.eqv;
%
%  Note:  dimensions are transposed for writing to file.
%
[npos, ndim] = size(crd);
[ncon, nnpe] = size(con);
%
dxfield = {['object "' fname '" class field']};
header  = {};
%
%-------------------- Positions.
%
objnam = 'nodal points'; % object name
objcls = 'array';        % object class
objtyp = 'float';        % object type
objrnk = 1;
objshp = ndim;
objitm = npos;
objfil = [fname '-positions.dat'];
%
objdsc = DXObjDesc(...
    objnam, objcls, objtyp, objrnk, objshp, objitm, objfil);
%
objatr = {};   % No attributes needed for positions
%
dxfield = {dxfield{:}, ...
           '   component "positions" "nodal points"' };
header = {header{:}, objdsc{:}, objatr{:}};
save(objfil, 'crd', '-ascii');
%
%-------------------- Connections.
%
objnam = 'connectivity';   % object name
objcls = 'array';          % object class
objtyp = 'int';            % object type
objrnk = 1;
objshp = nnpe;
objitm = ncon;
objfil = [fname '-connections.dat'];
%
objdsc = DXObjDesc(...
    objnam, objcls, objtyp, objrnk, objshp, objitm, objfil);
%
%   Attributes.
%
if (nnpe == 2) 
  etype = '"lines"';
elseif (nnpe == 3) 
  etype = '"triangles"';
elseif (nnpe == 4)
  etype = '"tetrahedra"';
elseif (nnpe == 8)
  etype = '"cubes"';
else
  error('Cannot determine element type.');
end
%
if (opts.ElementType) 
  etype = ['"', opts.ElementType, '"'];
end
%
objatr = {...
    ['attribute "element type" string ' etype],
    'attribute "ref" string "positions"'
         };   
%
dxfield = {dxfield{:}, ...
           '   component "connections" "connectivity"' };
header = {header{:}, ' ', objdsc{:}, objatr{:}};
WriteDataFile(objfil, {}, con, 'integer');
%
%-------------------- Nodal data.
%
lenargs = length(Ndata);
nobj    = lenargs/2;
Ndata   = reshape(Ndata, [2 nobj]);
%
for i=1:nobj
  %
  name = Ndata{1, i};
  valu = Ndata{2, i};
  if (~isempty(eqv))
    valu = ToAllNodes(valu, eqv); % checks dimension
  end
  %
  [m, n] = size(valu);
  %
  if (n == 1)              % transpose
    valu = valu';
    [m, n] = size(valu);
  end
  objnam = name;           % object name
  objcls = 'array';        % object class
  objtyp = 'float';        % object type
  %
  if (m == 1)
    objrnk = 0;
    objshp = 0;
  else
    objrnk = 1;
    objshp = m;
  end
  %
  objitm = n;
  objfil = [fname '-' name '.dat'];
  %
  objdsc = DXObjDesc(...
      objnam, objcls, objtyp, objrnk, objshp, objitm, objfil);
  %
  objatr = {'attribute "dep" string "positions"'};   
  %
  dxfield = {dxfield{:}, ...
             ['   component "' name '"  "' name '"']};
  header = {header{:}, ' ', objdsc{:}, objatr{:}};
  %
  valu = valu';
  save(objfil, 'valu', '-ascii');
  %
end
%
%-------------------- Elemental data.
%
lenargs = length(Edata);
nobj    = lenargs/2;
Edata   = reshape(Edata, [2 nobj]);
%
for i=1:nobj
  %
  name = Edata{1, i};
  valu = Edata{2, i};
  %
  [m, n] = size(valu);
  %
  if (n == 1)              % transpose
    valu = valu';
    [m, n] = size(valu);
  end
  if (size(valu, 2) ~= ncon)
    error('Element data has inconsistent size.')
  end
  objnam = name;           % object name
  objcls = 'array';        % object class
  objtyp = 'float';        % object type
  %
  if (m == 1)
    objrnk = 0;
    objshp = 0;
  else
    objrnk = 1;
    objshp = m;
  end
  %
  objitm = n;
  objfil = [fname '-' name '.dat'];
  %
  objdsc = DXObjDesc(...
      objnam, objcls, objtyp, objrnk, objshp, objitm, objfil);
  %
  objatr = {'attribute "dep" string "connections"'};   
  %
  dxfield = {dxfield{:}, ...
             ['   component "' name '"  "' name '"']};
  header = {header{:}, ' ', objdsc{:}, objatr{:}};
  %
  valu = valu';
  save(objfil, 'valu', '-ascii');
  %
end
%-------------------- Write header.
%
header = {header{:}, ' ', ...
          dxfield{:}, 'end'
          };
dxfile = [fname '.dx'];
WriteDataFile(dxfile, header, []);
%
%-------------------- Functions.
%
           
function objs = DXObjDesc(name, class, type, rank, shape, items, file)
% OBJSTRING - Construct object descriptor for DX data files.
%
%   objs is a cell array of strings
%   
onum = ['object '   '"' name '"'];
ocls = [' class '   class ];
otyp = [' type '    type];
ornk = [' rank '    int2str(rank)];
if (rank > 0)
  oshp = [' shape '   int2str(shape)];
else
  oshp = ' ';
end
oitm = [' items '   int2str(items) ' ascii'];
ofil = ['data file ' file];
%
objs = {[onum ocls otyp ornk oshp oitm], ofil};
