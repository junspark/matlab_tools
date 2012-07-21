function prod = MultMatArray(ma1, ma2)
% MultMatArray - Multiply arrays of matrices.
%   
%   USAGE:
%
%   prod = MultMatArray(ma1, ma2)
%
%   INPUT:
%
%   ma1 is m x n x l, 
%       an array of m x n matrices
%   ma2 is n x k x l, 
%       an array of n x k matrices
%  
%   OUTPUT:
%
%   prod is m x k x l, 
%        the array of matrix products
%
[m1 n1 l1] = size(ma1);
[n2 k2 l2] = size(ma2);
%
if (l1 ~= l2)
  error('mismatch on number of matrices')
end
%
if (n1 ~= n2)
  error('mismatch on internal matrix dimensions')
end
%%%%
%%%%  Make a sparse matrix, multipy and return to original format.
%%%%
%%%prod = MatArrayOfSparse(...
%%%    SparseOfMatArray(ma1)*SparseOfMatArray(ma2), ...
%%%    [m1 k2 l2]);
%%%%
%
%  Using the loop gave (slightly) faster results.
%
prod = zeros(m1, k2, l1);
for j=1:l1
  prod(:, :, j) = ma1(:, :, j)*ma2(:, :, j);
end
