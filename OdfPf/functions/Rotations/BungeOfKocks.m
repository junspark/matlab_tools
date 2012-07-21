function bunge = BungeOfKocks(kocks, units)
% BungeOfKocks - Bunge angles from Kocks angles.
%   
%   USAGE:
%
%   bunge = BungeOfKocks(kocks, units)
%
%   INPUT:
%
%   kocks is 3 x n,
%         the Kocks angles for the same orientations
%   units is a string,
%         either 'degrees' or 'radians'
%
%   OUTPUT:
%
%   bunge is 3 x n,
%         the Bunge angles for n orientations 
%
%   NOTES:
%
%   *  The angle units apply to both input and output.
%
if (nargin < 2)
  error('need two arguments: kocks, units')
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
bunge = kocks;
%
bunge(1, :) = kocks(1, :) + pi_over_2;
bunge(3, :) = pi_over_2 - kocks(3, :);
