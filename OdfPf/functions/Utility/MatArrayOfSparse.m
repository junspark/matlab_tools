function ma = MatArrayOfSparse(smat, mnl)
% MATARRAYOFSPARSE - Extract array of matrices from sparse format.
%   
%   smat is a j x k sparse matrix:
%           it is assumed to a block matrix with all blocks
%           the same shape
%   mnl  is an integer vector of length three:
%           it gives the shape of the block matrix `smat' (m x n x l)
%
%   ma is an m x n x l array of reals:
%         it is the array of m x n matrices derived from
%         the sparse input matrix
%
%   The input matrix must satisfy j=l*m and k=l*n.
%
xmax  = max(abs(smat(:))) + 1;
nzmat = SparseOfMatArray(repmat(xmax, [mnl(1) mnl(2) mnl(3)])) + smat;
ma    = reshape(nonzeros(nzmat), mnl) - xmax;

%old way% Significantly slower.
%old way% 
%old way% m = mnl(1); n = mnl(2); l = mnl(3);
%old way% for j=1:l
%old way%   ind_1 = (1:m) + (j-1)*m;
%old way%   ind_2 = (1:n) + (j-1)*n;
%old way%   ma(:, :, j) = full(smat(ind_1, ind_2));
%old way% end
