function rod = RodOfQuat(quat)
% RodOfQuat - Rodrigues parameterization from quaternion.
%   
%   USAGE:
%
%   rod = RodOfQuat(quat)
%
%   INPUT:
%
%   quat is 4 x n, 
%        an array whose columns are quaternion paramters; 
%        it is assumed that there are no binary rotations 
%        (through 180 degrees) represented in the input list
%
%   OUTPUT:
%
%  rod is 3 x n, 
%      an array whose columns form the Rodrigues parameterization 
%      of the same rotations as quat
% 
rod = quat(2:4, :)./repmat(quat(1,:), [3 1]);
