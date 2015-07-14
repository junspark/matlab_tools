function sv = SumValue(f, l2ip)
% SumValue - Find integral sum value of function.
%   
%   USAGE:
%
%   sv = SumValue(f, l2ip)
%
%   INPUT:
%
%   f    is an n-vector, 
%        an array of nodal point function values 
%   l2ip is n x n, 
%        the (L2) inner product matrix for the underlying mesh
%
%   OUTPUT:
%
%   sv is 1 x 1, 
%      the sum value of `f'---the integral of f over the region
%

sv  = sum(l2ip*f);
