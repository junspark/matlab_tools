function W = Skew(mat)

dim = size(mat, 1);

if dim ~= 3 | size(mat, 2) ~= dim
    error('your data is not 3-d, or it is not properly arranged!')
end

nmat = size(mat, 3);

matArray = SparseOfMatArray(mat);

skwArray = matArray - 0.5*(matArray + matArray');

W = MatArrayOfSparse(skwArray, [dim, dim, nmat]);
