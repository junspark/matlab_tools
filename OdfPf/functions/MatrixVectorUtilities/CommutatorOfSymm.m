function bigD = CommutatorOfSymm(d)
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
bigD = [...
  0                    ; -d(5, :)/r2          ;   d(6, :)/r2         ;...
	d(4, :)/r2           ;  0                   ;  -d(6, :)/r2         ;...
 -d(4, :)/r2           ;  d(5, :)/r2          ;   0                  ;...
 -(d(2, :)-d(3, :))/r2 ;  d(6, :)/2           ;  -d(5, :)/2          ;...
 -d(6, :)/2            ; (d(1, :)-d(3, :))/r2 ;   d(4, :)/2          ;...
  d(5, :)/2            ; -d(4, :)/2           ; -(d(1, :)-d(2, :))/r2;...
];

bigD = reshape(bigD, [3, 6, n]);
