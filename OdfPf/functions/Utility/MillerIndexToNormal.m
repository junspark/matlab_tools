function outvec = MillerIndexToNormal(invec, aspect)
% MILLERINDEXTONORMAL - Generate the normal(s) for a plane(s) given in 
%   by Miller indices in an triclinic basis {o1, o2, o3}.  The
%   basis for the output {c1, c2, c3} is chosen such that:
%
%       c1 || o1
%       c2 .  o2 > 0
%       c3 .  o3 > 0
%
% USAGE
%     outvec = MillerIndexToNormal(invec)
%
% INPUTS
%     1) invec: a 3 x n array of integers where n is the number of Miller
%         indices (hkl).
%
%   OUTPUTS
%     1) outvec: The 3 x n array whose columns are the unit vector
%         directions corresponding to the normal to the plane indices in 
%         the nth row/column of invec.
%   

if nargin == 1
  aspect = 1;
end

if length(aspect) == 1
  cob = diag([1, 1, 1/aspect]);
else
  cob = diag([1, 1/aspect(1), 1/aspect(2)]);    
end

outvec = UnitVector(cob*invec);

