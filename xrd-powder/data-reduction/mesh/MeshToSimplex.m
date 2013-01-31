function tm = MeshToSimplex(m, varargin)
% MeshToSimplex - Convert mesh to simplicial type.
%   
%   USAGE:
%
%    tm = MeshToSimplex(m)
%
%   INPUT:
%
%   m is a MeshStructure
%
%   OUTPUT:
%
%   tm is a MeshStructure
%      If `m' is already simplicial, then `tm' is the same.  If
%      `m' is rectangular, then each element of `m' is subdivided
%      in a regular way to produce `tm'.
%
%   NOTES:
%
%   * For cubes to tets, each cube is subdivided into six tets.
%
if m.etype.is_simplicial
  tm = m;
  warning('input element type is already simplicial');
end
%
intype = m.etype.name;
numelm = size(m.con, 2);
%
switch intype
  %
 case 'quads:4'
  newtype = 'triangles:3';
  econ = [
      1  2
      2  4
      3  3
      ];
  %
 case 'bricks:8'
  newtype = 'tets:4';
  econ = [
      1    2    4    4   4    4
      2    4    7    7   5    7
      3    3    3    5   2    6
      5    5    5    6   6    8
	  ];
  %
 otherwise
  error('conversion not available for type:  %s\n', intype)
end
%
esize = size(econ);
tm    = MeshStructure(m.crd, ...
		      reshape(m.con(econ, :), [esize(1),  esize(2)*numelm]), ...
		      [], 'ElementType', newtype);
%
return
