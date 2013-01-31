function MeshExportDX(m, fname, varargin)
% MeshExportDX - Export mesh data.
%   
%   USAGE:
%
%   MeshExportDX(m, fname, <options>)
%
%   INPUT:
%
%   m is a MeshStructure
%
%   fname is a string
%         The output filename.
%
%   OUTPUT:
%
%   *** NONE *** 
%
%   NOTES:
%
%  * Need to include option for passing integer data fields (e.g. grain number)
%
if isstr(fname)
  fid = fopen(fname, 'w');
else
  error('file input is neither integer nor string')
end
%
try
  DXEtype = DXElemType(m.etype.name);
catch
  error(lasterr);
end
MyName = mfilename;
%
%-------------------- Defaults and Keyword Arguments
%
optcell = {...
    'ElementData',  {{}}, ...
    'NodalData',    {{}}  ...
       };
%
options = OptArgs(optcell, varargin);
%
%-------------------- Main Execution
%
C = DXInterfaceConstants;
%
dxf = DXFile('Name', fname);
%
%  Write positions and connections.
%
if isfield(m, 'ndiv')
  counts = m.ndiv(end:-1:1) + 1;
  sizeA1 = size(m.transform.A, 1);
  deltas =  m.transform.A(:, end:-1:1)./repmat(counts - 1, [sizeA1, 1]);
  pos = DXGridPositions('Name',  'MeshPositions',   ...
			'Counts', counts', ...
			'Origin', m.transform.x0',   ...
			'Deltas', deltas);
  con = DXGridConnections('Name', 'MeshConnections', ...
			  'Counts', counts  ...
			  );
else
  pos = DXArray('Name', 'MeshPositions',   'Value', m.crd);
  con = DXArray('Name', 'MeshConnections', 'Value', m.con - 1, 'Type', 'int', ...
		'Attributes', [C.RefPositions, DXEtype]);
end
%
write(dxf, pos, con);
%
%  Write data arrays.
%
FieldHasData = (~isempty(options.ElementData) | ~isempty(options.NodalData));
%
if FieldHasData
  DataComponents = repmat(DXObjectRef(), [0 1]);
  DataNames      = {};
end
%
if ~isempty(options.ElementData)
  for i=1:2:length(options.ElementData)
    name = options.ElementData{i};
    val  = options.ElementData{i+1};
    dxa  = DXArray('Name', name,   'Value', val, 'Attributes', [C.DepConnections]);
    write(dxf, dxa);
    DataNames = {DataNames{:}, name};
    DataComponents = [DataComponents; DXObjectRef('Object', dxa,    'RefName', name)];
  end
end
%
if ~isempty(options.NodalData)
end
%
%  Write field.
%
FldComponents = [
    DXObjectRef('Object', pos,    'RefName', 'positions'),
    DXObjectRef('Object', con,    'RefName', 'connections')
    ];
FieldAtts = [];
if FieldHasData
  FldComponents = [FldComponents; DataComponents];
  %
  dxa = DXArray('Name', 'Field Names', 'Type', 'string', 'Value', DataNames);
  write(dxf, dxa);
  FieldAtts     = [DXAttribute('Name', 'Field Names', 'Kind', 'ObjectRef', ...
			 'Value', DXObjectRef('Object', dxa))];
end
%
DataField = DXField('Name', 'Generic Field', 'Components', FldComponents, ...
		    'Attributes', FieldAtts);
%
write(dxf, DataField);
%
close(dxf);
%
return
%
%grains = DXArray('Name', 'grains', 'Value', myspecs.grains, ...;
%		 'Attributes', [DEP_CON]);
%%
%FldComponents = [
%    DXObjectRef('Object', pos,    'RefName', 'positions'),
%    DXObjectRef('Object', con,    'RefName', 'connections')
%    DXObjectRef('Object', grains, 'RefName', 'Grain Number')
%    ];
%
%%
%write(dxf, grains);
%
function t = DXElemType(intype)
% DXElemType - Return DX Element Type attribute
%
switch intype
 case 'lines:2'
  t = 'lines';
 case 'triangles:3'
  t = 'triangles';
 case 'tets:4'
  t = 'tetrahedra';
 case 'quads:4'
  t = 'quads';
 case 'bricks:8'
  t = 'bricks';
 otherwise
  error('no DX output available for type %s\n', intype)
end
%
t = DXAttribute('Name', 'element type', 'Value', t);
%
return
