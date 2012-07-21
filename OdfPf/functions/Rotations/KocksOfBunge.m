function kocks = KocksOfBunge(bunge, units)
% KocksOfBunge - Kocks angles from Bunge angles.
%   
%   USAGE:
%
%   kocks = KocksOfBunge(bunge, units)
%
%   INPUT:
%
%   bunge is 3 x n,
%         the Bunge angles for n orientations 
%   units is a string,
%       either 'degrees' or 'radians'
%
%   OUTPUT:
%
%   kocks is 3 x n,
%         the Kocks angles for the same orientations
%
%   NOTES:
%
%   *  The angle units apply to both input and output.
%
if (nargin < 2)
  error('need two arguments: bunge, units')
end
%
if (strcmp(units, 'degrees'))
  %
  indeg = 1;
  %
elseif (strcmp(units, 'radians'))
  %
  indeg = 0;
  %
else
  error('angle units need to be specified:  ''degrees'' or ''radians''')
end
%
if (indeg)
  pi_over_2 = 90;
else
  pi_over_2 = pi/2;
end
%
kocks = bunge;
%
kocks(1, :) = bunge(1, :) - pi_over_2;
kocks(3, :) = pi_over_2 - bunge(3, :);
