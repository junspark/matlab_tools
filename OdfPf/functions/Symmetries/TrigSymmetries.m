function [varargout] = TrigSymmetries()
% HEXSYMMETRIES - Generate quaternions for the hexagonal symmetry group as
%   well as the parameters describing the facets of the fundamental region
%   in Rodrigues' space.
%
%   Usage:
%   
%       1.            trig_sym = TrigSymmetries
%       2. [trig_sym, FRFaces] = TrigSymmetries
%
%   Inputs:
%
%       none
%
%   Outputs:
%
%       1. trig_sym: The 4 x n array of unit quaternions representing the 
%           trigonal symmetry group where n is the number of symmetries
%           (6).
%
%       2. FRFaces*: The 4 x m array whose columns represent 
%           [phi_symmetry, [face normal]]' and m is the number of
%           unique faces (4).
%   
  AngleAxis = [...
      0.0      0    0    0    1; % identity
      pi*2/3   0    0    0    1; % threefold about c 
      pi*4/3   0    0    0    1;
      pi       1   -2    1    0; % twofold about 2(-1)(-1)0
      pi       2   -1   -1    0;
      pi       1    1   -2    0;
	      ]';
  %
  AngleAxis(2:4, :) = MillerBravaisToUnit(AngleAxis(2:5, :), 1);
  %
  if nargout == 2
    temp = [...
        pi*2/3   0    0    0    1;
	pi       1   -2    1    0; % twofold about 2(-1)(-1)0
	pi       2   -1   -1    0;
	pi       1    1   -2    0;
	   ]';
    %
    FRFaces(1, :) = temp(1, :);
    FRFaces(2:4, :) = MillerBravaisToUnit(temp(2:5, :), 1);
    varargout{2} = FRFaces;
  end
  %
  Angle = AngleAxis(1,:);
  Axis  = AngleAxis(2:4,:);
  %
  %  Axis does not need to be normalized in call to QuatOfAngleAxis.
  %
  varargout{1} = QuatOfAngleAxis(Angle, Axis);
  