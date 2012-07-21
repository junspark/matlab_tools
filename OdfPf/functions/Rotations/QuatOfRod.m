function quat = QuatOfRod(rod)
% QuatOfRod - Quaternion from Rodrigues vectors.
%   
%   USAGE:
%
%   quat = QuatOfRod(rod)
%
%   INPUT:
%
%   rod  is 3 x n, 
%        an array whose columns are Rodrigues parameters
%
%   OUTPUT:
%
%   quat is 4 x n, 
%        an array whose columns are the corresponding unit
%        quaternion parameters; the first component of 
%        `quat' is nonnegative
%
cphiby2   = cos(atan(sqrt(dot(rod,rod,1))));
%
quat = [cphiby2 ; repmat(cphiby2, [3 1]) .* rod];
