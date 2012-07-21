function mat = Decomp2Matrix(vec)
% usage:
%  mat = Decomp2Matrix(vec)
%  
%   vec is 6 x n
%   mat is 3 x 3 x n (symm) 
num_vec = size(vec, 2);
%
mat = zeros(3, 3, num_vec);
for i = 1:num_vec
    mat(:, :, i) = diag(vec(1:3, i)) + 1/sqrt(2)*(diag([vec(6, i), vec(4, i)], 1) + diag(vec(5, i), 2));
    mat(:, :, i) = mat(:, :, i) + triu(mat(:, :, i), 1)';
end
