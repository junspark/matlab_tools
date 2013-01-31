function m10 = MeshT4ToT10(m4)
%  MeshT4ToT10 - Convert 4 to 10 noded tets.
%
%  *** Recently changed order of refpts to match tets:10 definition
%      This alters FemEvps output.
%   
%   USAGE:
%
%   m10 = MeshT4ToT10(m4)
%
%   INPUT:
%
%   m4 is a MeshStructure
%      with element type 'tets:4'
%
%   OUTPUT:
%   
%   m10 is a MeshStructure
%       with element type 'tets:10'
%
%   NOTES:
%
%   * Does not handle equivalences 
%   * Output mesh is unsorted 
%   * Could be extended to handle triangles:3 --> triangles:6 and lines:2 --> lines:3
%
if ~strcmp(m4.etype.name, 'tets:4')
  error('input mesh is not marked as four-node tets')
end
%
tol   = 1.0e-7;
numel = size(m4.con, 2);
%
refpts = [
    1.0  0.0  0.0  0.0
    0.5  0.5  0.0  0.0  % 1,2 (midpoint)
    0.0  1.0  0.0  0.0
    0.5  0.0  0.5  0.0  % 1,3
    0.0  0.5  0.5  0.0  % 2,3
    0.0  0.0  1.0  0.0
    0.5  0.0  0.0  0.5  % 1,4
    0.0  0.5  0.0  0.5  % 2,4
    0.0  0.0  0.5  0.5  % 3,4
    0.0  0.0  0.0  1.0
	 ]';
%
%  Create new nodal points including duplicates.
%
newcrds = SpreadRefPts(m4, refpts);
%
%  Create new connectivity.
%
newcon = reshape(1:(10*numel), [10 numel]);
%
%  Find unique positions.
%
[newcrds, ord, iord] = UniqueVectors(newcrds, tol);
%
%  Adjust the connectivity for the new node numbers.
%
newcon = iord(newcon);
%
m10 = MeshStructure(newcrds, newcon, [], 'ElementType', 'tets:10');
%
%old refpts% refpts = [1.0  0.0  0.0  0.0     % 1  (new number)
%old refpts% 	  0.5  0.5  0.0  0.0  % 1,2   % 2
%old refpts% 	  0.0  1.0  0.0  0.0          % 3
%old refpts% 	  0.0  0.5  0.5  0.0  % 2,3   % 5
%old refpts% 	  0.0  0.0  1.0  0.0          % 6
%old refpts% 	  0.5  0.0  0.5  0.0  % 1,3   % 4
%old refpts% 	  0.5  0.0  0.0  0.5  % 1,4   % 7
%old refpts% 	  0.0  0.5  0.0  0.5  % 2,4   % 8
%old refpts% 	  0.0  0.0  0.5  0.5  % 3,4   % 9
%old refpts% 	  0.0  0.0  0.0  1.0          % 10
%old refpts% 	 ]';
