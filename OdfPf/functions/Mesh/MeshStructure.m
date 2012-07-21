function mesh = MeshStructure(crd, con, eqv, varargin)
% MeshStructure - Create mesh structure from mesh data.
%   
%   USAGE:
%
%   mesh = MeshStructure
%   mesh = MeshStructure(crd, con)
%   mesh = MeshStructure(crd, con, eqv)
%   mesh = MeshStructure(crd, con, eqv, <options>)
%
%   INPUT:
%
%   crd is e x n
%       the array of nodal point locations
%   con is d x m  (integer)
%       the mesh connectivity
%   eqv is 2 x k, (integer, optional) 
%       the equivalence array
%
%   <options> is a sequence of parameter, value pairs
%             Available options are listed below with default values
%             shown in brackets.
%
%   'ElementType'    string
%                    See ElementTypeStruct for available types.
%   'Symmetries'     4 x n
%                    an array of quaternions for crystal symmetries
%
%   OUTPUT:
%
%   mesh is a MeshStructure, 
%        the basic MeshStructure consists of three fields,
%        the nodal point coordinates (.crd), the connectivity
%        (.con) and the nodal point equivalence array (.eqv).
%        
%   NOTES:
% 
%   *  With no arguments, this function returns and empty mesh
%      structure.  With only two arguments, it sets the equivalence 
%      array to be empty.
%
%
MyName = mfilename;
%
%-------------------- Defaults and Keyword Arguments
%
optcell = {...
    'ElementType',  '', ...
    'Symmetries',   []  ...
       };
%
options = OptArgs(optcell, varargin);
%
%-------------------- Main Code
%
if (nargin == 0)
  crd = [];
  con = [];
  eqv = [];
elseif (nargin == 2)
  eqv = [];
end
%
mesh = struct('crd', crd, 'con', con, 'eqv', eqv);
%
%  Process optional arguments.
%
if nargin > 3
  %
  if options.ElementType
    mesh.etype = ElementTypeStruct(options.ElementType);
  end
  %
  if options.Symmetries
    mesh.symmetries = options.Symmetries;
  end
  %
end
