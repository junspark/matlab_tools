function rmat = RMatOfQuat(quat)
% RMatOfQuat - Convert quaternions to rotation matrices.
%   
%   USAGE:
%
%   rmat = RMatOfQuat(quat)
%
%   INPUT:
%
%   quat is 4 x n, 
%        an array of quaternion parameters
%
%   OUTPUT:
%
%   rmat is 3 x 3 x n, 
%        the corresponding array of rotation matrices
%
%   NOTES:
%
%   *  This is not optimized, but still does okay
%      (about 6,700/sec on intel-linux ~2GHz)
%
n    = size(quat, 2);
rmat = zeros(3, 3, n);
%
zerotol = 1.0e-7;  % sqrt(eps) due to acos
for i=1:n
  theta = 2*acos(quat(1, i));
  if (theta > zerotol)                % find axis
    w = (theta/sin(theta/2))*quat(2:4, i);
  else
    w = [0 0 0]';
  end
  %
  wskew = [0 -w(3) w(2);  w(3) 0 -w(1); -w(2) w(1) 0];
  %
  rmat(:, :, i) = expm(wskew);
  %
end
