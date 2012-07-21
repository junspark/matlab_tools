function q = ToFundamentalRegionQ(quat, qsym)
% ToFundamentalRegionQ - To quaternion fundamental region.
%   
%   USAGE:
%
%   q = ToFundamentalRegionQ(quat, qsym)
%
%   INPUT:
%
%   quat is 4 x n, 
%        an array of n quaternions
%   qsym is 4 x m, 
%        an array of m quaternions representing the symmetry group
%
%   OUTPUT:
%
%   q is 4 x n, the array of quaternions lying in the
%               fundamental region for the symmetry group 
%               in question
%
%   NOTES:  
%
%   *  This routine is very memory intensive since it 
%      applies all symmetries to each input quaternion.
%
n = size(quat, 2); % number of points
m = size(qsym, 2); % number of symmetries
%
%  Apply all symmetries to each member of quat.
%  
qeqv = QuatProd(...
    reshape(repmat(quat, m, 1), 4, m*n), ...
    repmat(qsym, 1, n));
%
%  Calculate max cosine and select corresponding element.
%
[qmax, imax] = max(abs(reshape(qeqv(1, :), m, n)), [], 1);
indices = (0:n-1)*m + imax;
q = qeqv(:, indices);
