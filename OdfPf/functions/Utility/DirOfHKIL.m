function dir = DirOfHKIL(hkil, cbya)
% DIROFHKIL - Compute Cartesian direction from Miller-Bravais indices.
%   
%   dir = DirOfHKIL(hkil)
%   dir = DirOfHKIL(hkil, cbya)
%  
%   hkil is a 4 x m real array:
%           it is an array of m sets of Miller-Bravais indices
%   cbya is a real scalar:  (optional) 
%           it is the c/a ratio for the hexagonal system; 
%           default value is one
%
%   dir is a 3 x m real array:
%          it is the list of unit vectors corresponding
%          to the specified Miller-Bravais indices
%    
if (nargin < 2)
  cbya = 1.0;
end
%
%  Check that first three indices sum to zero.
%
check = [1 1 1] * hkil(1:3, :);
if (max(abs(check)) > 1.0e-6)
  error('First three indices do not sum to zero.')
end
%
angs = (0:2)*(2*pi/3);
c    = cos(angs);
s    = sin(angs);
z    = [0 0 0];
%
vecs = [[c 0]; [s 0]; [z cbya]];
%
dir = UnitVector(vecs*hkil);
