function qp = QuatProd(q2, q1)
% QuatProd - Product of two unit quaternions.
%   
%   USAGE:
%
%    qp = QuatProd(q2, q1)
%
%   INPUT:
%
%    q2, q1 are 4 x n, 
%           arrays whose columns are quaternion parameters
%
%   OUTPUT:
%
%    qp is 4 x n, 
%       the array whose columns are the quaternion parameters of 
%       the product; the first component of qp is nonnegative
%
%    NOTES:
%
%    *  If R(q) is the rotation corresponding to the
%       quaternion parameters q, then 
%       
%       R(qp) = R(q2) R(q1)
%
a = q2(1, :); a3 = repmat(a, [3, 1]);
b = q1(1, :); b3 = repmat(b, [3, 1]);
%
avec = q2(2:4, :);
bvec = q1(2:4, :);
%
qp = [...
    a.*b - dot(avec, bvec);                ...
    a3.*bvec + b3.*avec + cross(avec, bvec)...
    ];
%
q1negative = (qp(1,:) < 0 );
qp(:, q1negative) = -1*qp(:, q1negative);
