function S = Symm(mat)

dim = size(mat, 1);

if dim ~= 3 | size(mat, 2) ~= dim
    error('your data is not 3-d, or it is not properly arranged!')
end

nmat = size(mat, 3);

matArray = SparseOfMatArray(mat);

symArray = 0.5*(matArray + matArray');

S = MatArrayOfSparse(symArray, [dim, dim, nmat]);
