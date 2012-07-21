function values = EvalMeshFunc(mesh, fun, els, ecrds)
% EvalMeshFunc - Evaluate mesh function.
%   
%   USAGE:
%
%   values = EvalMeshFunc(mesh, fun, els, ecrds)
%
%   INPUT:
%
%   mesh  is a MeshStructure
%   fun   is an n-vector, 
%         function values on the set of independent nodal points
%   els   is an m-vector of integers, 
%         it gives the list of elements containing the points at
%         which to evaluate the function
%   ecrds is k x m, 
%         it is the list of barycentric coordinates of the 
%         points at which to evaluate the function
%
%   OUTPUT:
%
%   values is an n-vector,
%          the values of the function at the specified points
%
fun = fun(:)';  % make row vector
eqv = mesh.eqv;
if (~isempty(eqv))
  fun = ToAllNodes(fun, eqv);
end
%
con  = mesh.con;
fcon = fun(con(:, els));
%
values = dot(fcon, ecrds, 1);

