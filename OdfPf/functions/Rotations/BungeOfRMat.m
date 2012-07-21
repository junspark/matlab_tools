function bunge = BungeOfRMat(rmat, units)
% BungeOfRMat - Bunge angles from rotation matrices.
%   
%   USAGE:
%
%   bunge = BungeOfRMat(rmat, units)
%
%   INPUT:
%
%   rmat  is 3 x 3 x n, 
%         an array of rotation matrices
%   units is a string,
%         either 'degrees' or 'radians' specifying the output
%         angle units
%
%   OUTPUT:
%
%   bunge is 3 x n,
%         the array of Euler angles using Bunge convention
%
if (nargin < 2)
  error('need second argument, units:  ''degrees'' or ''radians''')
end
%
if (strcmp(units, 'degrees'))
  %
  indeg = 1;
  %
elseif (strcmp(units, 'radians'))
  indeg = 0;
else
  error('angle units need to be specified:  ''degrees'' or ''radians''')
end
%
warning off MATLAB:divideByZero; % allow divide by zero and fix it later
%
n  = size(rmat, 3);
%
c2 = rmat(3, 3, :);
c2 = min(c2(:), 1.0);
c2 = max(c2', -1.0);
%
myeps         = sqrt(eps);  %  for arc-cosine
near_pole     = (abs(c2) > 1-myeps);
%not_near_pole = ~near_pole;
%
%  Not near pole.
%
s2 = sqrt(1 - c2.*c2);
%
c1 = -squeeze(rmat(2, 3, :))' ./s2; % squeeze makes column vector here
s1 =  squeeze(rmat(1, 3, :))' ./s2;
c3 =  squeeze(rmat(3, 2, :))' ./s2;
s3 =  squeeze(rmat(3, 1, :))' ./s2;
%
%  Near pole.
%
c1(near_pole) = squeeze(rmat(1, 1, near_pole))';
s1(near_pole) = squeeze(rmat(2, 1, near_pole))';
c3(near_pole) = 1;
s3(near_pole) = 0;
%
bunge = [atan2(s1, c1); acos(c2); atan2(s3, c3)];
bunge(bunge < 0) = bunge(bunge < 0) + 2*pi;
%
if (indeg)
  bunge = bunge * (180/pi);
end
