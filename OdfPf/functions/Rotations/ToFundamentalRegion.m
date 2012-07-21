function rod = ToFundamentalRegion(quat, qsym)
% ToFundamentalRegion - Put rotation in fundamental region.
%   
%   USAGE:
%
%   rod = ToFundamentalRegion(quat, qsym)
%
%   INPUT:
%
%   quat is 4 x n, 
%        an array of n quaternions
%   qsym is 4 x m, 
%        an array of m quaternions representing the symmetry group
%
%   OUTPUT:
%
%   rod is 3 x n, 
%       the array of Rodrigues vectors lying in the fundamental 
%       region for the symmetry group in question
%
%   NOTES:  
%
%   *  This routine is very memory intensive since it 
%      applies all symmetries to each input quaternion.
%
q   = ToFundamentalRegionQ(quat, qsym);
rod = RodOfQuat(q);
