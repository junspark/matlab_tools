function bigD = CommutatorOfSymm2(d)
% bigD = CommutatorOfSymm(d)
%
% d    -- vector of a symmetric matrix (6 x n)
% bigD -- Skew operator in the matrix-vector notation
%         to perform [D]*[x] - [x]*[D] as [bigD]*{x} where [D] is
%         the symmetric matrix of {d} and [x] and {x} are the matrix and
%         vector representations of a symmetric tensor.
% 

if size(d, 1) ~= 6
    error('your data is not 6-d, or it is not arranged column-wise!')
end
r2 = sqrt(2);

n = size(d, 2);
bigD = [0 ; r2*d(4, :) ; -r2*d(4, :) ; r2*(d(3, :)-d(2, :)) ; -d(6, :) ; d(5, :) ;...
        -r2*d(5, :) ; 0 ; r2*d(5, :) ; d(6, :) ; r2*(d(1, :)-d(3, :)) ; -d(4, :) ;...
        r2*d(6, :) ; -r2*d(6, :) ; 0 ; -d(5, :) ; d(4, :) ; r2*(d(2, :)-d(1, :)) ;...
];

bigD = reshape(bigD, [6, 3, n]);
