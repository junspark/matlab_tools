function osym = OrtSymmetries()
% OrtSymmetries - Orthorhombic symmetry group.
%
%   USAGE:
%
%   osym = OrtSymmetries
%
%   INPUT:  none
%
%   OUTPUT:
%
%   osym is 4 x 4, 
%        the quaternions for the symmetry group
%   
AngleAxis = [...
    0.0      1    1    1; % identity
    pi       1    0    0;
    pi       0    1    0;
    pi       0    0    1;
	    ]';
%
Angle = AngleAxis(1,:);
Axis  = AngleAxis(2:4,:);
%
%  Axis does not need to be normalized in call to QuatOfAngleAxis.
%
osym = QuatOfAngleAxis(Angle, Axis);
