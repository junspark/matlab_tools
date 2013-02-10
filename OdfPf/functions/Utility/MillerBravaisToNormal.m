function outvec = MillerBravaisToNormal(invec, aspect)
% MILLERBRAVAISTONORMAL - Generate the normal(s) for a plane(s) given in
%   the Miller-Bravais convention for the hexagonal basis {a1, a2, a3, c}.  The
%   basis for the output {o1, o2, o3} is chosen such that:
%
%       o1 || a1
%       o3 || c
%       o2 = o3 ^ o1
%
%   Usage:
%
%       outvec = MillerBravaisToUnit(invec, aspect)
%
%   Inputs:
%
%       1. invec: Either a 4 x n or n x 4 array of integers
%           where n is the number of Miller-Bravais directions [uvtw]
%           s.t. u, v, t, w are a reduced set of integers with u + v + t = 0.
%       2. aspect (optional): is the 1 x 1 or 2 x 1 scalar array of the aspect
%           ratio(s).
%
%   Outputs:
%
%       1. outvec: The 3 x n array whose columns are the unit vector
%           directions corresponding to the normal to the plane indices in
%           the nth row/column of invec.
%

if nargin == 1
    aspect = 1;
end

size_vec = size(invec);

if size_vec(1) == 4 & size_vec(2) ~= 4
    cols = size_vec(2);
elseif size_vec(2) == 4 & size_vec(1) ~= 4
    cols = size_vec(1);
    invec = invec';
elseif size_vec(1) == 4 & size_vec(2) == 4
    if sum(sum(invec(1:3, :))) == 0
        cols = size_vec(2);
    elseif sum(sum(invec(:, 1:3))) == 0
        cols = size_vec(1);
        invec = invec';
    else
        error('There is a problem with your input: At least one of the inputs violates u + v + t = 0');
        return
    end
else
    error('There is a problem with your input: The input list is dimensionally incorrect')
    return
end

outvec = UnitVector(...
    [invec(1, :);...
    (2*invec(2, :) + invec(1, :))/sqrt(3);...
    invec(4, :)*1/aspect;...
    ]);
