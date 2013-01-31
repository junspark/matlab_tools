function [m, activenodes] = MeshClean(m0)
% MeshClean - Eliminate unused nodes.
%   
%   USAGE:
%
%   m = MeshClean(m0)
%   [m, activenodes] = MeshClean(m0)
%
%   INPUT:
%
%   m0 is a MeshStructure
%      It may have nodes that are not referenced by the connectivity
%
%   OUTPUT:
%
%   m           is a MeshStructure
%     Any unreferenced nodes are eliminated and the mesh is renumbered.
%   activenodes is an integer array
%               It is the list node numbers in the original mesh which
%               are active.
%
%   NOTES:
%
%   *  It would be useful to add a 'parent_nodes' list to the output
%
m = m0;
activenodes = unique(m.con(:));
numnodes    = length(activenodes);
newnum      = sparse(ones(1, numnodes), activenodes, 1:numnodes);
m.crd       = m.crd(:, activenodes);
m.con       = full(newnum(m.con));
%
%  For one element case, matlab indexing is not as expected.
%  (always returns a row vector)
%
s0 = size(m0.con);
s1 = size(m.con);
if any(s1 - s0)
  m.con = reshape(m.con, s0);
end
%
return
