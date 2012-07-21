function bigW = CommutatorOfSkew(w)
% bigW = CommutatorOfSkew(w)
%
% w    -- axial vector of skew matrix (3 x n)
% bigW -- Skew operator in the matrix-vector notation
%         to perform [W]*[x] - [x]*[W] as [bigW]*{x} where [W] is
%         the skew matrix of {w} and [x] and {x} are the matrix and
%         vector representations of a symmetric tensor.
% 

if size(w, 1) ~= 3
    error('your data is not 3-d, or it is not arranged column-wise!')
end
r2 = sqrt(2);

n = size(w, 2);
bigW = [zeros(4, n);-r2*w(2, :);r2*w(3, :);...
        zeros(3, n);r2*w(1, :);zeros(1, n);-r2*w(3, :);...
        zeros(3, n);-r2*w(1, :);r2*w(2, :);zeros(1, n);...
        zeros(1, n);-r2*w(1, :);r2*w(1, :);zeros(1, n);-w(3, :);w(2, :);...
        r2*w(2, :);zeros(1, n);-r2*w(2, :);w(3, :);zeros(1, n);-w(1, :);...
        -r2*w(3, :);r2*w(3, :);zeros(1, n);-w(2, :);w(1, :);zeros(1, n);
];

bigW = reshape(bigW, [6, 6, n]);

