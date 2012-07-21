function quat = QuatRotate(qa, q0)
% QUATROTATE - Apply a change of basis matrix to quaternion array.
%   
%   quat = QuatRotate(qa, q0)
%
%   qa is 4 x n, an array of n quaternions
%   q0 is 4 x 1, a single quaternion
%
%   quat is 4 x n, the array of rotated quaternions
%
%   for each quaternion Q in qa, the result is
%   Q0 * Q * Q0^T
%
n    = size(qa, 2);
q0n  = repmat(q0, [1 n]); % n copies of q0
quat = QuatProd(q0n, qa); % multiply on left
%
%  Now apply inverse.
%
q0n(1, :) = -q0n(1, :);   % n copies of q0 inverse
quat = QuatProd(quat, q0n); % multiply on right
