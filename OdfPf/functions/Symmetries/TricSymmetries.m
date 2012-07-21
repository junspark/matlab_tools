function sym = TricSymmetries()
% TRICSYMMETRIES - Generate quaternions for triclinic symmetry group.
%
%   hsym = HexSymmetries
%
%   hsym is a 4 x 12 real array:
%           it is the hexagonal symmetry group represented
%           as quaternions
%   
angleAxis = [...
        0,   1,   0,   0,...
]';
% 
Angle = angleAxis(1, :);
Axis  = angleAxis(2:4, :);
%
%  Axis does not need to be normalized in call to QuatOfAngleAxis.
%
sym = QuatOfAngleAxis(Angle, Axis);
