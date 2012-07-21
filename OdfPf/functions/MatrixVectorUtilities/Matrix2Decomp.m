function vec = Matrix2Decomp(mat)
% usage:
%  vec = Matrix2Decomp(mat)
%  
%   vec is 6 x n
%   mat is 3 x 3 x n (symm) 
num_mat = size(mat, 3);
%
vec = zeros(6, num_mat);
for i = 1:num_mat
    sym_test = triu(mat(:, :, i))' - tril(mat(:, :, i));
%     if norm(sym_test, inf) < eps
%         disp(['input is not symmetric', num2str(i)])
%         return
%     end
    vec(1:3, i) = diag(mat(:, :, i));
    vec(4, i) = sqrt(2)*mat(2, 3, i);
    vec(5, i) = sqrt(2)*mat(1, 3, i);
    vec(6, i) = sqrt(2)*mat(1, 2, i);
end
