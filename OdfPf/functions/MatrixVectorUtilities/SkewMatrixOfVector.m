function W = SkewMatrixOfVector(w)
% SKEWMATRIXOFVECTOR - Find the skew matrix(matrices) associated with the given axial
% vector(s)
%
% USAGE
%     W = SkewMatrixOfVector(w)
%
% INPUTS
%     1) w is 3 x n, an array of n horizontally concatenated axial vectors
%        in R^3.
%
% OUTPUTS
%     1) W is 3 x 3 x n, the array of n skew matrices associated with each
%        axial vector.
%
% SEE ALSO
%     VectorOfSkewMatrix

if size(w, 1) ~= 3
  error('your data is not 3-d, or it is not arranged column-wise!')
end

dim3 = size(w, 2);

W = [...
     zeros(1, dim3);
            w(3, :);
           -w(2, :);
           -w(3, :);
     zeros(1, dim3);
            w(1, :);
            w(2, :);
           -w(1, :);
     zeros(1, dim3);
];

W = reshape(W, [3, 3, dim3]);

%% Or, if you want it to run slowly... %%

% % % W = zeros(3, 3, dim3);
% % % for i = 1:dim3
% % %     W(:, :, i) = W(:, :, i) + diag(-1.0*[w(3, i), w(1, i)], 1);
% % %     W(:, :, i) = W(:, :, i) + diag(w(2, i), 2);
% % %     W(:, :, i) = W(:, :, i) - triu(W(:, :, i))';
% % % end
