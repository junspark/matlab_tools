function minv = InvMatArray(marray)
% InvMatArray - Inverses for an array of matrices.
%   
%   USAGE:
%
%   minv = InvMatArray(marray)
%
%   INPUT:
%
%   marray is n x n x l, 
%          an array of n x n matrices
%
%   OUTPUT:
%
%   minv is n x n x l, 
%        the inverse of each matrix
%
[n n1 l] = size(marray);
%
if (n ~= n1)
  error('matrices not square');
end
%
minv = zeros(n, n, l);
for i = 1:l
  minv(:, :, i) = inv(marray(:, :, i));
end
%
