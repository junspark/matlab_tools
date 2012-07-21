function hsym = HexSymmetries()
% HexSymmetries - Quaternions for hexagonal symmetry group.
%
%   USAGE:
%
%   hsym = HexSymmetries
%
%   INPUT:  none
%
%   OUTPUT:
%
%   hsym is 4 x 12,
%        it is the hexagonal symmetry group represented
%        as quaternions
%   
p3  = pi/3;
p6  = pi/6;
ci  = cos(p6*(0:5));
si  = sin(p6*(0:5));
z6  = zeros(1, 6);
w6  = ones (1, 6);
pi6 = repmat(pi, [1 6]);
%
sixfold = [ [p3*(0:5)]; z6; z6; w6];
twofold = [pi6; ci; si; z6];
%
AngleAxis = [sixfold twofold];
%
Angle = AngleAxis(1,:);
Axis  = AngleAxis(2:4,:);
%
%  Axis does not need to be normalized in call to QuatOfAngleAxis.
%
hsym = QuatOfAngleAxis(Angle, Axis);
