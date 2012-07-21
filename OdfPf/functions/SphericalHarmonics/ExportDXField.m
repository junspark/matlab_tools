function fref = ExportDXField(fname, mesh, Ndata, Edata, varargin)
% ExportDX1 - Write a DX field in DX native format.
%   
%   USAGE:
%
%   header = ExportDX(fname, mesh)
%   header = ExportDX(fname, mesh, Ndata)
%   header = ExportDX(fname, mesh, Ndata, Edata)
%   header = ExportDX(fname, mesh, Ndata, Edata, <options>)
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
%                     per element, 'tetrahedra' is chosen over 
%                     'quads'; possible values are:
%                     'lines'|'triangles'|'quads'|'tetrahedra'|'cubes'
%   'FieldName'       string specifying name of DX field 
%                     defaults to 'Field'
%
%   OUTPUT:
%
%   fref is a DXObjectRef
%           referring to the output field of this routine
%
%   NOTES:
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
Ndata = {Ndata{:}};
Edata = {Edata{:}};
try 
  Ndata = reshape(Ndata, [2 length(Ndata)/2]);
  Edata = reshape(Edata, [2 length(Edata)/2]);
catch
  error('Odd number of entries in either Ndata or Edata');
end

%
%-------------------- Defaults and Keyword Arguments
%
optcell = {...
    'ElementType',  '',       ...
    'FieldName',    'Field'   ...
	  };
%
opts = OptArgs(optcell, varargin);
%
etype = opts.ElementType;
if isempty(etype);
  etype = getElementType(size(mesh.con, 1)); % determine from number of nodes per element
end
%
%-------------------- Standard Attributes
%
REF_POS  = DXAttribute('Name', 'ref', 'Value', 'positions');
DEP_POS  = DXAttribute('Name', 'dep', 'Value', 'positions');
DEP_CON  = DXAttribute('Name', 'dep', 'Value', 'connections');
%
%-------------------- Open DXFile
%
fname = deblank(fname);
if ~strcmp(fname((end-2):end), '.dx' )
  fname = [fname, '.dx'];
end
dxf = DXFile('Name', fname);
%
%-------------------- Write the mesh.
%
ATT_ETYPE = DXAttribute('Name', 'element type', 'Value', etype);
crd = DXArray('Name', 'Nodal Points', ...
	      'Value', mesh.crd);
con = DXArray('Name', 'Elements', ...
	      'Value', mesh.con - 1, ...
	      'Type', 'int', ...
	      'Attributes', [REF_POS, ATT_ETYPE]);
%
MeshComponents = [
    DXObjectRef('Object', crd, 'RefName', 'positions'),
    DXObjectRef('Object', con, 'RefName', 'connections') 
		 ]';
%
dxf = write(dxf, crd, con);
clear crd con
%
%-------------------- Nodal Data
%
nobj    = size(Ndata, 2);
NDFlds  = cell(1, nobj);
NDComps = repmat(DXObjectRef(), [1 nobj]);
if nobj
  for i=1:nobj
    refName = Ndata{1, i};
    dtype   = getDataType(Ndata{2, i});
    %
    %  Apply nodal equivalences if necessary.
    %
    if (isempty(mesh.eqv))
      values = Ndata{2, i};
    else
      values = ToAllNodes(Ndata{2, i}, mesh.eqv);
    end
    %
    dxa = DXArray('Name', refName, ...
		  'Value', values, ...
		  'Type', dtype, ...
		  'Attributes', [DEP_POS]);
    dxf = write(dxf, dxa);
    NDFlds{i} = refName;
    NDComps(i) = DXObjectRef('Object', dxa);
  end
end
%
%-------------------- Elemental Data
%
nobj    = size(Edata, 2);
EDFlds  = cell(1, nobj);
EDComps = repmat(DXObjectRef(), [1 nobj]);
if nobj
  for i=1:nobj
    refName = Edata{1, i};
    dtype   = getDataType(Edata{2, i});
    dxa = DXArray('Name', refName, ...
		  'Value', Edata{2, i}, ...
		  'Type', dtype, ...
		  'Attributes', [DEP_CON]);
    dxf = write(dxf, dxa);
    EDFlds{i}  = refName;
    EDComps(i) = DXObjectRef('Object', dxa);
  end
end
%
%-------------------- Write Field
%
DataFields = {NDFlds{:}, EDFlds{:}};
if length(DataFields)
  DataFNames = DXArray('Name', 'Data Components', 'Type', 'string', 'Value', DataFields);
  DCAtt = DXAttribute('Name', 'Data Components', ...
		      'Kind', 'ObjectRef', ...
		      'Value', DXObjectRef('Object', DataFNames));
  write(dxf, DataFNames);
  FldAtts = [DCAtt];
else
  FldAtts = [];
end
dxfld = DXField('Name', opts.FieldName, ...
		'Components', [MeshComponents, NDComps, EDComps],...
		'Attributes', FldAtts);
dxf  = write(dxf, dxfld);
fref = DXObjectRef('Object', dxfld, 'File', fname);
%
%-------------------- Helper Functions.
%
function etype = getElementType(nnpe)
% getElementType - Element type from number of nodes per element.
%   
%
switch nnpe
 case 2
  etype = 'lines';
 case 3
  etype = 'triangles';
 case 4
  etype = 'tetrahedra';
 case 8
  etype = 'cubes';
 otherwise
  error('Cannot determine element type.');
end

return

function dtype = getDataType(x)
% getDataType - determine DX datatype from matlab datatype
%
%  Should take place in DXInterface module in @DXArray
%
switch class(x)
 case 'single'
  dtype = 'float';
 case 'double'
  dtype = 'double';
 case 'char'
  dtype = 'string';
 case {'int8', 'unint8', 'int16', 'uint16', 'int32', 'unint32', 'int64', 'uint64'}
  dtype = 'int';
 otherwise
  error(['Array data type not matched:  ', class(x)])
end
