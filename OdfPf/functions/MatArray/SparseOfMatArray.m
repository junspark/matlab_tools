function smat = SparseOfMatArray(marray)
% SparseOfMatArray - Sparse matrix from array of matrices.
%   
%   USAGE:
%
%   smat = SparseOfMatArray(marray)
%
%   INPUT:
%
%   marray is m x n x l, 
%          an array of m x n matrices
%
%   OUTPUT:
%
%   smat is l*m x l*n , (sparse)
%        the sparse matrix form of the block matrix `marray'
%
[m n l] = size(marray);
%
mn   = m*n;
jmax = l*n;
imax = l*m;
ntot = l*m*n;
%
sij = reshape(marray, [1 ntot]);
j   = reshape(repmat( (1:jmax), [m 1]), [1 ntot]);
i   = reshape(repmat( (1:m), [1 jmax]), [1 ntot]) + ...
      reshape(repmat( (m*(0:(l-1))),   [mn 1]), [1 ntot]);
%
smat = sparse(i, j, sij, imax, jmax);
