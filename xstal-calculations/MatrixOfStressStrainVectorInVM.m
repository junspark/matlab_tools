function t = MatrixOfStressStrainVectorInVM(v, varargin)
% MatrixOfStressStrainVectorInVM - Converts vectorized stress (strain) to a
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
%   OUTPUT:
%
%   t
%       stress (strain) matrix.
%
%   NOTE:
%       This method of vectorizing and associated rotation works as long 
%       as the shear components do not contribute to the normal components
%       and vice versa.

% default options
optcell = {...
    'Order', '11-22-33-23-13-12', ...
    };

% update option
opts    = OptArgs(optcell, varargin);

v   = v(:);
if length(v) == 6
    s2  = sqrt(2);
    if strcmp(opts.Order, '11-22-33-23-13-12')
        t   = [ ...
            v(1)    v(6)/s2 v(5)/s2; ...
            v(6)/s2    v(2) v(4)/s2; ...
            v(5)/s2 v(4)/s2 v(3); ...
            ];
    elseif strcmp(opts.Order, '11-22-33-12-13-23')
        t   = [ ...
            v(1)    v(4)/s2 v(5)/s2; ...
            v(4)/s2    v(2) v(6)/s2; ...
            v(5)/s2 v(6)/s2 v(3); ...
            ];
    else
        disp('order not supported');
        return
    end
else
    disp('input vector size incorrect');
end