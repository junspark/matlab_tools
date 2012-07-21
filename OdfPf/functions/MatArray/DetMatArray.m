function determ = DetMatArray(marray)
% DetMatArray - Evaluate determinants for an array of matrices.
%   
%   USAGE:
%
%   determ = DetMatArray(marray)
%
%   INPUT:
%
%   marray is n x n x l, 
%          an array of l n x n matrices
%
%   OUTPUT:
%
%   determ is 1 x l, 
%          the determinant of each matrix
%
[n n1 l] = size(marray);
%
if (n ~= n1)
  error('matrices not square');
end
%
determ = zeros(1, l);
for i = 1:l
  determ(i) = det(marray(:, :, i));
end
%
%  Below is code for using the sparse LU decomposition,
%  but it was not faster for larger problems.
%
%disp('using lu')
%tic
%[L, U, P] = lu(SparseOfMatArray(marray));
%udiag  = reshape(diag(U), [n l]);
%dlu    = prod(udiag, 1);
%toc
%%
