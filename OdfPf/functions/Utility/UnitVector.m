function uvec = UnitVector(vec, ipmat)
% UnitVector - Normalize an array of vectors.
%   
%   USAGE:
%
%   uvec = UnitVector(vec)
%   uvec = UnitVector(vec, ipmat)
%
%   INPUT:
%
%   vec   is m x n, 
%         an array of n nonzero vectors of dimension m
%   ipmat is m x m, (optional)
%         this is a (SPD) matrix which defines the inner product
%         on the vectors by the rule:  
%            norm(v)^2 = v' * ipmat * v
%         
%         If `ipmat' is not specified, the usual Euclidean 
%         inner product is used.
%
%   OUTPUT:
%
%   uvec is m x n,
%        the array of unit vectors derived from `vec'
%
m = size(vec, 1);
%
if (nargin > 1) 
  nrm2 = dot(vec, ipmat*vec, 1)
else
  nrm2 = dot(vec, vec, 1);
end
%
nrm  = repmat(sqrt(nrm2), [m 1]);
uvec = vec./nrm;
