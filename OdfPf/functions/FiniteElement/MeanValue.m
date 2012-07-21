function mv = MeanValue(f, l2ip)
% MeanValue - Find integral mean value of function.
%   
%   USAGE:
%
%   mv = MeanValue(f, l2ip)
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
%   mv is 1 x 1, 
%      the mean value of `f'---the integral of f divided by
%      the measure of the region
%
n     = length(f);
nones = ones(1, n);
%
iwts = full(sum(l2ip));
mv = iwts*f(:)./sum(iwts);
