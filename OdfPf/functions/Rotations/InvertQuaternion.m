function invq = InvertQuaternion(q)
% INVERTQUATERNION - Return inverse quaternions for a list of unit quaternions
%
% USAGE
%     1) invq = InvertQuaternion(q)
%
% INPUTS
%     1) q is 4 x n, a list of n horizontally concatenated unit quaternions.
%
% OUTPUTS
%     1) invq is 4 x n, a list of n horizontally concatenated unit
%         quaternions representing the inverse of each point in q; i.e. the
%         quaternion product (invq * q) generates the identity
%
% SEE ALSO
%     QuatProd, QuatOfAngleAxis, QuatOfRod, QuatRotate
%
% NOTES
%     1) Yes, semi-silly utility routine
[check, n] = size(q);
%
if check ~= 4
    error('Either input is not a list of quaternions or is not horizontally concatenated');
end
%
% Simply put...
% I choose to leave the angle positive (right-handed)
invq = [q(1, :);-q(2:4, :)];
