function rmat = RMatOfKocks(kocks, units)
% RMATOFKOCKS - Rotation matrix from Kocks angles.
%   
%   rmat = RMatOfKocks(kocks, units)
%
%   kocks is 3 x n,
%         the array of Kocks parameters (in radians)
%   units is a string,
%         either 'degrees' or 'radians'
%       
%   rmat is 3 x 3 x n,
%	the corresponding rotation matrices
%
if (nargin < 2)
  error('need second argument, units:  ''degrees'' or ''radians''')
end
%
if (strcmp(units, 'degrees'))
  %
  indeg = 1;
  kocks = kocks*(pi/180);
  %
elseif (strcmp(units, 'radians'))
  indeg = 0;
else
  error('angle units need to be specified:  ''degrees'' or ''radians''')
end

numRot = size(kocks, 2);

cosPSI = cos(kocks(1, :));
cosTHT = cos(kocks(2, :));
cosPHI = cos(kocks(3, :));

sinPSI = sin(kocks(1, :));
sinTHT = sin(kocks(2, :));
sinPHI = sin(kocks(3, :));

% Rotation of CRYSTAL -> SAMPLE
% [Kocks, et. al., "Texture and Anisotropy", Cambridge, UK 1998]
rmat = [
    -sinPHI.*sinPSI - cosPHI.*cosPSI.*cosTHT;
     sinPHI.*cosPSI - cosPHI.*sinPSI.*cosTHT;
     cosPHI.*sinTHT;
     cosPHI.*sinPSI - sinPHI.*cosPSI.*cosTHT;
    -cosPHI.*cosPSI - sinPHI.*sinPSI.*cosTHT;
     sinPHI.*sinTHT;
     cosPSI.*sinTHT;
     sinPSI.*sinTHT;
     cosTHT;
       ];

rmat = reshape(rmat, [3, 3, numRot]);
