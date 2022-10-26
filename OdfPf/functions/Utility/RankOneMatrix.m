function r1mat = RankOneMatrix(vec1, vec2)
% RankOneMatrix - Create rank one matrices (dyadics) from vectors.
%   
%   USAGE:
%
%   r1mat = RankOneMatrix(vec1)
%   r1mat = RankOneMatrix(vec1, vec2)
%
%   INPUT:
%
%   vec1 is m1 x n, 
%        an array of n m1-vectors 
%   vec2 is m2 x n, (optional) 
%        an array of n m2-vectors
%
%   OUTPUT:
%
%   r1mat is m1 x m2 x n, 
%         an array of rank one matrices formed as c1*c2' 
%         from columns c1 and c2
%
%   With one argument, the second vector is taken to
%   the same as the first.
%
%   NOTES:
%
%   *  This routine can be replaced by MultMatArray.
%
if (nargin == 1)
  vec2 = vec1;
end
%
[m1, n1] = size(vec1);
[m2, n2] = size(vec2);
%
if (n1 ~= n2)
  error('Number of points differ in arguments.')
end
%
m1m2 = m1 * m2;
%
r1mat = zeros(m1m2, n1);
%
range = 1:m1;
%
for i=1:m2
  r1mat(range, :) = vec1 .* repmat(vec2(i, :), [m1 1]);
  range = range + m1;
end
%
r1mat = reshape(r1mat, [m1 m2 n1]);
