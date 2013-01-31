function eqv = DuplicateVectors(vec, tol)
% DUPLICATEVECTORS - Find vectors in an array that are equivalent to within
% a specified tolerance
%
%   USAGE:
%
%       eqv = DuplicateVectors(vec, <tol>)
%
%   INPUT:
%
%       1) vec is n x m, a double array of m horizontally concatenated
%                        n-dimensional vectors.
%      *2) tol is 1 x 1, a scalar tolerance.  If not specified, the default
%                        tolerance is 1e-14.
%
%   OUTPUT:
%
%       1) eqv is 1 x p, a cell array of p equivalence relationships.
%
%   NOTES:
%
%       Each equivalence relationship is a 1 x q vector of indices that
%       represent the locations of duplicate columns/entries in the array
%       vec.  For example:
%
%             | 1     2     2     2     1     2     7 |
%       vec = |                                       |
%             | 2     3     5     3     2     3     3 |
%
%       eqv = {[1x2 double]    [1x3 double]}, where
%
%       eqv{1} = [1     5]
%       eqv{2} = [2     4     6]

eqv = [];

if nargin == 1
    tol = 1e-14;
end

vlen = size(vec, 2);
orid  = [1:vlen];

torid = orid;
tvec  = vec;

ii = 1;
while vlen > 1
    dupl = repmat(tvec(:, 1), [1, vlen]);

    diff = sum(abs(tvec - dupl), 1);

    match = find(abs(diff(:, 2:end)) <= tol) + 1;

    kick = [1, match];

    if length(kick) > 1;
        eqv{ii} = torid(kick);
        ii = ii + 1;
    end

    cmask       = ones(1, vlen);
    cmask(kick) = 0;
    cmask       = logical(cmask);

    tvec  = tvec(:, cmask);

    torid = torid(cmask);

    vlen = size(tvec, 2);
end

