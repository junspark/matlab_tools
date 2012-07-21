function new = RotConvert(old, from, to)
% ROTCONVERT - Convert rotation representation.
%   
%   STATUS:  in development
%
%   USAGE:
%
%   new = RotConvert(old, from, to)
%
%   INPUT:
%
%   old  is m1 x n
%        a list of n orientations in some format
%   from is a string
%        It names the format of the input rotations ("old")
%   to   is a string
%        It specifies the output format.
%
%   OUTPUT:
%
%   new is m2 x n
%       It is the same rotations represented by "old", but using the "to" format
%
%   NOTES:
%
%   * expects degrees for Euler angles
%
r   = ConvertFrom(old, from);
new = ConvertTo(r, to);
return
%
%  ---------------------------------- Helper functions
%
function n = ConvertFrom(old, from)
% CONVERTFROM - return rotation matrices
%   
%
switch lower(from)
 case 'rmat'
  n = old;
 case 'quat'
  n = RMatOfQuat(old);
 case 'bunge'
  n = RMatOfBunge(old, 'degrees');
 case 'kocks'
  n = RMatOfBunge(BungeOfKocks(old, 'degrees'), 'degrees');
 case 'rod'
  n = RMatOfQuat(QuatOfRod(old));
 otherwise
  error('"from" argument not recognized')
end
return

%
function n = ConvertTo(r, to)
% CONVERTFROM - return rotation matrices
%   
%
switch lower(to)
 case 'rmat'
  n = r;
 case 'quat'
  n = QuatOfRMat(r);
 case 'bunge'
  n = BungeOfRMat(r, 'degrees');
 case 'kocks'
  n = KocksOfBunge(BungeOfRMat(r, 'degrees'), 'degrees');
 case 'rod'
  n = RodOfQuat(QuatOfRMat(r));
 otherwise
  error('"to" argument not recognized')
  
end
return

