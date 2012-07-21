function csym = CubSymmetries()
% CubSymmetries - Return quaternions for cubic symmetry group.
%
%   USAGE:
%
%   csym = CubSymmetries
%
%   INPUT:  none
%
%   OUTPUT:
%
%   csym is 4 x 24, 
%        quaternions for the cubic symmetry group
%
AngleAxis = [...
    0.0      1    1    1; % identity
    pi*0.5   1    0    0; % fourfold about x1
    pi       1    0    0;
    pi*1.5   1    0    0;
    pi*0.5   0    1    0; % fourfold about x2
    pi       0    1    0;
    pi*1.5   0    1    0;
    pi*0.5   0    0    1; % fourfold about x3
    pi       0    0    1;
    pi*1.5   0    0    1;
    pi*2/3   1    1    1; % threefold about 111
    pi*4/3   1    1    1;
    pi*2/3  -1    1    1; % threefold about 111
    pi*4/3  -1    1    1;
    pi*2/3   1   -1    1; % threefold about 111
    pi*4/3   1   -1    1;
    pi*2/3  -1   -1    1; % threefold about 111
    pi*4/3  -1   -1    1;
    pi       1    1    0; % twofold about 110
    pi      -1    1    0;
    pi       1    0    1;
    pi       1    0   -1;
    pi       0    1    1;
    pi       0    1   -1;
	    ]';
%
Angle = AngleAxis(1,:);
Axis  = AngleAxis(2:4,:);
%
%  Axis does not need to be normalized; it is done
%  in call to QuatOfAngleAxis.
%
csym = QuatOfAngleAxis(Angle, Axis);
