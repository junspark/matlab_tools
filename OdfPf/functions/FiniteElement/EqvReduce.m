function rmat = EqvReduce(mat, eqv)
% EqvReduce - Reduce matrix columnwise by equivalences.
%   
%   USAGE:
%
%   rmat = EqvReduce(mat, eqv)
%
%   INPUT:
%
%   mat is m x n, (usually sparse) a matrix, the columns of which
%                 have underlying nodal point equivalences
%   eqv is 2 x k, the equivalence array
%
%   OUTPUT:
%
%   rmat is m x (n-k), the new matrix for the reduced (condensed) 
%                 set of nodes formed by adding equivalent columns 
%                 to the master column
%
%   NOTES:
%
%   *  The columns of the unreduced nodes are added to the 
%      columns of the master node.
%
%   *  This routine only reduces along the column dimension.  To
%      reduce along the rows, apply to the transpose.  To reduce
%      along both dimensions, call it twice.
%
%   *  The equivalence array can be empty, in which case nothing
%      is done, and the original matrix is returned.
%
if (isempty(eqv))
  rmat = mat;
  return
end
%
n = size(mat, 2);
k = size(eqv, 2);
j = n - k;
%
for i=1:k;
  iu = eqv(1, i); 
  ir = eqv(2, i);
  %
  mat(:, ir) = mat(:, ir) + mat(:, iu);
end
%
rmat = mat(:, 1:j);
