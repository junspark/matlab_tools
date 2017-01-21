function t = VectorOfStressStrainMatrixInVM(v)
% MatrixOfStressStrainInVM - Converts vectorized stress (strain) to a
% matrix
%
%   USAGE:
%
%   t = MatrixOfStressStrainInVM(v)
%
%   INPUT:
%
%   v
%       modified Voight-Mandel stress (strain) vector in [11, 22, 33, 12,
%       13, 23] order. Shears are expected to have factor of sqrt(2)
%       multiplied to them already.
%
%   t
%       stress (strain) matrix.

if (size(v,1) == 3) && (size(v,2) == 3)
    s2  = sqrt(2);
    t   = [ ...
        v(1,1); ...
        v(2,2); ...
        v(3,3); ...
        s2*v(1,2); ...
        s2*v(1,3); ...
        s2*v(2,3); ...
        ];
else
    disp('input matrix size incorrect');
end