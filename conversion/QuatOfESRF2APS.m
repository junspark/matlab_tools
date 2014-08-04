function quat = QuatOfESRF2APS
% QuatOfESRF2APS - Quaternion that converts ESRF coordinates to APS
% coordinates
%   
%   USAGE:
%
%   quat = QuatOfESRF2APS
%
%   INPUT:
%
%   None
%
%   OUTPUT:
%
%   quat is 4 x n,
%        the quaternion representation of the change of basis matrix that
%        takes ESRF coordinate system to APS coordinate system
%
%   Note:
%
%   ESRF coordinate system
% 
%       X - along x-ray beam
%       Y - OB
%       Z - UP
%
%   APS coordinate system
%       X - OB
%       Y - UP
%       Z - along the x-ray beam

RMat    = [ ...
    0 1 0; ...
    0 0 1; ...
    1 0 0; ...
    ];
quat    = QuatOfRMat(RMat);
