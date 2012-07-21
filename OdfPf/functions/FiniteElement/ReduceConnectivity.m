function conr = ReduceConnectivity(mesh)
% ReduceConnectivity - Apply equivalences to connectivity.
%   
%   USAGE:
%
%   conr = ReduceConnectivity(mesh)
%
%   INPUT:
%
%   mesh is a MeshStructure
%
%   OUTPUT:
%
%   conr is m x n, 
%        the connectivity array with values all
%        referring to the set of independent nodes; if the
%        mesh.eqv array is empty, conr is the same as 
%        mesh.con
%
con = mesh.con;
eqv = mesh.eqv;
%
if(isempty(eqv))
  conr = con;
  return
end
%
nnp  = max(con(:));
neq  = size(eqv, 2);
nred = nnp - neq;
%
[e1, i1, j1] = unique(eqv(1, :));
ord1 = [1:nred, nred+i1];
con = ord1(con);
eqv2 = eqv(2, i1);
ord2 = [1:nred eqv(2, :)];
conr = ord2(con);
%
