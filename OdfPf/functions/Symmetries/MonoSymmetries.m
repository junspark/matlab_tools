function qsym = MonoSymmetries()
% MONOSYMMETRIES - Quaternions for Monoclinic symmetry group.
%
%   osym = OrtSymmetries
%
%   osym is 4 x n, where n is the number of symmetries (4)
%   
AngleAxis = [...
    0.0      1    1    1; % identity
    pi       1    0    0;
	    ]';
%
Angle = AngleAxis(1,:);
Axis  = AngleAxis(2:4,:);
%
%  Axis does not need to be normalized in call to QuatOfAngleAxis.
%
qsym = QuatOfAngleAxis(Angle, Axis);
