function RMat = RMatOfExpMap(ExpMap)
% RMATOFEXPMAP - Calculate the rotation matrix of one or more
% exponential map parameterization(s).
%   
%   USAGE:
%           RMat = RMatOfExpMap(ExpMap)
%           
%   INPUTS:
%           ExpMap is 3 x n, the exponential map parameter(s) for n
%           orientations arranged column-wise.
%   
%   OUTPUTS:
%           RMat is 3 x 3 x n, the rotation matrix corresponding to
%           each of the n inputs concatenated in the 3rd dimension.
%   
%   NOTES:
%           *) First rev 04/02/2004, uses a loop over the
%           closed-form expression for expm(skew(w))...this could
%           be replaced by a vectorized form to improve efficiency
%           if necessary.
%   
%   SEE ALSO:  RMatOfQuat, RMatOfBunge, RMatOfKocks, RMatOfRoe
  
  numObj = size(ExpMap, 2);
  
  % Angles
  mag = NormVecArray(ExpMap);
  
  % Skew matrices of axial vectors
  SkewMat = SkewMatrixOfVector(ExpMap);
  
  % Find zero angles to apply limits
  zeroIndex = find(mag < sqrt(eps));
  
  mag(zeroIndex) = 1;
  
  C1 = sin(mag)./mag;
  C1(zeroIndex) = 1;
  
  C2 = (1 - cos(mag))./mag.^2;
  C2(zeroIndex) = 1;
  
  % Form rotation matrices from closed-form expression for
  % expm(Skew(w))
  %
  % *) May repalce the loop with a vectorized form, although the
  % typical number of arguments may not justify it...
  Iden = eye(3);
  RMat = zeros(3, 3, numObj);
  for i = 1:numObj
    RMat(:, :, i) = Iden ...
	+ C1(i)*SkewMat(:, :, i) ...
	+ C2(i)*SkewMat(:, :, i)*SkewMat(:, :, i);
  end
  