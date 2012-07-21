function quat = QuatOfRMat(rmat)
% QuatOfRMat - Quaternion from rotation matrix
%   
%   USAGE:
%
%   quat = QuatOfRMat(rmat)
%
%   INPUT:
%
%   rmat is 3 x 3 x n,
%        an array of rotation matrices
%
%   OUTPUT:
%
%   quat is 4 x n,
%        the quaternion representation of `rmat'

% 
%  Find angle of rotation.
%
ca = 0.5*(rmat(1, 1, :) + rmat(2, 2, :) + rmat(3, 3, :) - 1);
ca = min(ca, +1);
ca = max(ca, -1);
angle = squeeze(acos(ca))';
%
%  Three cases for the angle:  
%  
%  *   near zero -- matrix is effectively the identity
%  *   near pi   -- binary rotation; need to find axis
%  *   neither   -- general case; can use skew part
%
tol = 1.0e-4;
anear0 = (angle < tol);
nnear0 = length(anear0);
angle(anear0) = 0;

raxis = [rmat(3, 2, :) - rmat(2, 3, :);
	rmat(1, 3, :) - rmat(3, 1, :);
	rmat(2, 1, :) - rmat(1, 2, :)];
raxis = squeeze(raxis);
raxis(:, anear0) = 1;
%
special = angle > pi - tol;
nspec   = sum(special);
if (nspec > 0)
  %disp(['special: ', num2str(nspec)]);
  sangle = repmat(pi, [1, nspec]);
  tmp = rmat(:, :, special) + repmat(eye(3), [1, 1, nspec]);
  tmpr = reshape(tmp, [3, 3*nspec]);
  tmpnrm = reshape(dot(tmpr, tmpr), [3, nspec]);
  [junk, ind] = max(tmpnrm);
  ind = ind + (0:3:(3*nspec-1));
  saxis = squeeze(tmpr(:, ind));
  raxis(:, special) = saxis;
end
%debug.angle = angle;
%debug.raxis = raxis;
quat = QuatOfAngleAxis(angle, raxis);

