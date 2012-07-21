function [varargout] = OrtSymmetries()
% ORTSYMMETRIES - Quaternions for orthorhombic symmetry group.
%
%   osym = OrtSymmetries
%
%   osym is 4 x n, where n is the number of symmetries (4)
%   
AngleAxis = [...
    0.0      1    1    1; % identity
    pi       1    0    0;
    pi       0    1    0;
    pi       0    0    1;
	    ]';
%
if nargout == 2
    temp = [...
        pi       0    0    1;
        pi       0    1    0;
        pi       1    0    0;
    ]';
    %
    FRFaces(1, :) = temp(1, :);
    FRFaces(2:4, :) = UnitVector(temp(2:4, :));
    varargout{2} = FRFaces;
end
%
Angle = AngleAxis(1,:);
Axis  = AngleAxis(2:4,:);
%
%  Axis does not need to be normalized in call to QuatOfAngleAxis.
%
varargout{1} = QuatOfAngleAxis(Angle, Axis);
