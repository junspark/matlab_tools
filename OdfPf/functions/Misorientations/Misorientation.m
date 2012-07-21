function [angle, mis] = Misorientation(q1, q2, sym)
% Misorientation - Return misorientation data for quaternions.
%
%   USAGE:
%
%   angle = Misorientation(q1, q2, sym)
%   [angle, mis] = Misorientation(q1, q2, sym)
%
%   INPUT:
%
%   q1 is 4 x n1, 
%      is either a single quaternion or a list of n quaternions
%   q2 is 4 x n,  
%      a list of quaternions
% 
%   OUTPUT:
%
%   angle is 1 x n, 
%         the list of misorientation angles between q2 and q1
%   mis   is 4 x n, (optional) 
%         is a list of misorientations in the fundamental region 
%         (there are many equivalent choices)
%
%   NOTES:
%
%   *  The misorientation is the linear tranformation which
%      takes the crystal basis given by q1 to that given by
%      q2.  The matrix of this transformation is the same
%      in either crystal basis, and that is what is returned
%      (as a quaternion).  The result is inverse(q1) * q2.
%      In the sample reference frame, the result would be
%      q2 * inverse(q1).  With symmetries, the result is put
%      in the fundamental region, but not into the Mackenzie cell.
%
[f1, n1] = size(q1);
[f2, n2] = size(q2);
%
if (n1 == 1)
  q1 = repmat(q1, [1 n2]);
end
%
q1i = [-q1(1, :); q1(2:4, :)];
mis = ToFundamentalRegionQ(QuatProd(q1i, q2), sym);
%
angle = 2*acos(min(1, mis(1, :)));
%
