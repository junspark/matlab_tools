function csym = CubSymmetries()
% CUBSYMMETRIES - Generate quaternions for cubic symmetry group.
%
%   csym = CubSymmetries
%
%   csym is 4 x n, where n is the number of symmetries (24)
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
%  Axis does not need to be normalized in call to QuatOfAngleAxis.
%
csym = QuatOfAngleAxis(Angle, Axis);
