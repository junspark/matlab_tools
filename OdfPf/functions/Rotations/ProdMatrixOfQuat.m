function qmats = ProdMatrixOfQuat(quats, flag)
% PRODMATRIXOFQUAT - Form 4 x 4 matrices to perform the quaternion
% product
%
% USAGE
%    qmats = ProdMatrixOfQuat(quats, flag)
%
% INPUTS
%    1) quats is 4 x n, an array of n quaternions horizonatlly concatenated
%    2) flag is a string, etiher 'left' or 'right' denoting the sense of
%       the multiplication, e.g.
%
%                    / ProdMatrixOfQuat(h, 'right') * q
%        q * h  --> <
%                    \ ProdMatrixOfQuat(q, 'left') * h      
%
% OUTPUTS
%    1) qmats is 4 x 4 x n, the left or right quaternion product operator
%
nq = size(quats, 2);

if strcmp(flag, 'right')
    temp = [...
        quats(1, :);
        quats(2, :);
        quats(3, :);
        quats(4, :);
        %
       -quats(2, :);
        quats(1, :);
       -quats(4, :);
        quats(3, :);
        %
       -quats(3, :);
        quats(4, :);
        quats(1, :);
       -quats(2, :);
        %
       -quats(4, :);
       -quats(3, :);
        quats(2, :);
        quats(1, :);
        ];
elseif strcmp(flag, 'left')
    temp = [...
        quats(1, :);
        quats(2, :);
        quats(3, :);
        quats(4, :);
        %
       -quats(2, :);
        quats(1, :);
        quats(4, :);
       -quats(3, :);
        %
       -quats(3, :);
       -quats(4, :);
        quats(1, :);
        quats(2, :);
        %
       -quats(4, :);
        quats(3, :);
       -quats(2, :);
        quats(1, :);
        ];
end

temp = reshape(temp, [4, 4, nq]);

for i = 1:nq
    qmats{i} = temp(:, :, i);
end
