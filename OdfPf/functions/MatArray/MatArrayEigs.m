function [evals, evecs] = MatArrayEigs(ma)
% MatArrayEigs - Eigenvalues of matrix array.
%   
%   USAGE:
%
%   [evals, evecs] = MatArrayEigs(ma)
%
%   INPUT:
%
%   OUTPUT:
%
%   NOTES:
%
%   * matlab returns eigenvalues from smallest to largest
%
sma = size(ma);
%
evals = zeros(sma(1), sma(3));
evecs = zeros(sma);
%
for i=1:sma(3)
  [evecs(:, :, i), D] = eig(ma(:, :, i));
  evals(:, i) = diag(D);
end
