function m4 = MeshT10ToT4(m10)
%  MeshT10ToT4 - Convert 10-node tets to 4-node tets
%
%   
%   USAGE:
%
%   m4 = MeshT10ToT4(m10)
%
%   INPUT:
%
%   m10 is a MeshStructure
%      with element type 'tets:10'
%
%   OUTPUT:
%   
%   m4 is a MeshStructure
%       with element type 'tets:4'
%
%   NOTES:
%
%   * Each 10-node tet produces 8 4-node tets
%   * Does not handle equivalences 
%   * Output mesh is unsorted 
%
if ~strcmp(m10.etype.name, 'tets:10')
  error('input mesh is not marked as 10-node tets')
end
%
tol   = 1.0e-7;
numel = size(m10.con, 2);
%
T4CON = [
    8     2     4     2     9     7     7     4
    2     9     9     4     7     2     4     7
    5     8     5     9     8     9     2     2
    3     5     6     5    10     8     1     9
    ];
%
%  Create new connectivity.
%
newcon = reshape(m10.con(T4CON, :), [4, 8*numel]);
%
m4 = MeshStructure(m10.crd, newcon, [], 'ElementType', 'tets:4');
%
