function w = VectorOfSkewMatrix(W)

if size(W, 1) ~= 3
    error('your data is not 3-d, or it is not arranged column-wise!')
end

dim2 = size(W, 3);
w = zeros(3, dim2);
for i = 1:dim2
    w(:, i) = [-W(2, 3, i), W(1, 3, i), -W(1, 2, i)]';
end
