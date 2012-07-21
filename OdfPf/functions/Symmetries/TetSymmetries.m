function varargout = TetSymmetries()
% TETSYMMETRIES - Generate quaternions for cubic symmetry group.
%
%   sym = CubSymmetries
%
%   sym is 4 x n, where n is the number of symmetries (24)
%   
  AngleAxis = [...
      0.0      1    1    1; % identity
      pi*0.5   0    0    1; % fourfold about x3 (001)
      pi       0    0    1;
      pi*1.5   0    0    1;
      pi       1    0    0; % twofold about x1 (100)
      pi       0    1    0; % twofold about x2 (010)
      pi       1    1    0; % twofold about 110
      pi      -1    1    0;
  ]';
  %
  Angle = AngleAxis(1,:);
  Axis  = AngleAxis(2:4,:);
  %
  if nargout == 2
    FRFaces = [...
	pi*0.5   0    0    1;
	pi       1    0    0;
	pi       0    1    0;
	pi       1    1    0; % twofold about 110
	pi      -1    1    0;
	      ]';
    %
    temp = FRFaces(2:4, :);
    temp = UnitVector(temp);
    FRFaces(2:4, :) = temp;
    %
    varargout{2} = FRFaces;
  end
  %
  Angle = AngleAxis(1,:);
  Axis  = AngleAxis(2:4,:);
  %
  %  Axis does not need to be normalized in call to QuatOfAngleAxis.
  %
  varargout{1} = QuatOfAngleAxis(Angle, Axis);
  