function p = ProjectionVector(n)
% ProjectionVector - Generate a projection vector that projects a
% vectorized 2nd order symmetric tensor when 2nd order symmetric tensor is
% vectorized as the following:
%         
%   [A] = [A_11, A_12, A_13;
%          A_12, A_22, A_23;
%          A_13, A_23, A_33]
%       |
%       |
%       V
%   {A} = [A_11, A_22, A_33, sqrt(2)*A_12, sqrt(2)*A_13, sqrt(2)*A_23]
%           
%   where the operation n*A*n' (in tensor notation) is obtained by the
%   matrix-vector product {g}*{A}.
%
%   Usage:
%   
%       g = ProjectionVector(n)
%
%   Input:
%
%       n   projection direction - unit vector (3 x 1)
%
%   Output:
%
%       p   (6 x 1) matrix that projects {A} in {n} 
% 
%   Note:
%   strain and stress are vectorized as the following:
%   [strain]    = [e11, e12, e13;
%                  e12, e22, e23;
%                  e13, e23, e33]
%
%   {strain}    = {e11
%                  e22
%                  e33
%                  sqrt(2)*e12
%                  sqrt(2)*e13
%                  sqrt(2)*e23}
% 
%   [stress]    = [s11, s12, s13;
%                  s12, s22, s23;
%                  s13, s23, s33]
%
%   {stress}    = {s11
%                  s22
%                  s33
%                  sqrt(2)*s12
%                  sqrt(2)*s13
%                  sqrt(2)*s23}
%
p   = [...
    n(1)^2 ...
    n(2)^2 ...
    n(3)^2 ...
    sqrt(2)*n(1)*n(2) ...
    sqrt(2)*n(1)*n(3) ...
    sqrt(2)*n(2)*n(3)];