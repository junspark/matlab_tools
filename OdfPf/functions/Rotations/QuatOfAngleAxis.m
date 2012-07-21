function quat = QuatOfAngleAxis(angle, raxis)
% QuatOfAngleAxis - Quaternion of angle/axis pair.
%
%   USAGE:
%
%   quat = QuatOfAngleAxis(angle, rotaxis)
%
%   INPUT:
%
%   angle is an n-vector, 
%         the list of rotation angles
%   raxis is 3 x n, 
%         the list of rotation axes, which need not
%         be normalized (e.g. [1 1 1]'), but must be nonzero
%
%   OUTPUT:
%
%   quat is 4 x n, 
%        the quaternion representations of the given
%        rotations.  The first component of quat is nonnegative.
%   
halfangle = 0.5*angle(:)';
cphiby2   = cos(halfangle);
sphiby2   = sin(halfangle);
%
rescale = sphiby2 ./sqrt(dot(raxis, raxis, 1));
%
quat = [cphiby2; repmat(rescale, [3 1]) .* raxis ] ;
%
q1negative = (quat(1,:) < 0);
quat(:, q1negative) = -1*quat(:, q1negative);
