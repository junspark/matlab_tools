function dN_dxi = SimplexSFunDer(dim)
% SIMPLEXSFUNDER -
%          | N1,1   N1,2   ....   N1,n |
%          | N2,1   N2,2   ....     .  |
% dN_dxi = |   .      .    .        .  |
%          |   .      .      .      .  |
%          | Nm,1     .        .   Nm,n|

dN_dxi = eye(dim);
dN_dxi = cat(1, dN_dxi, -1.0*ones(1, dim));
