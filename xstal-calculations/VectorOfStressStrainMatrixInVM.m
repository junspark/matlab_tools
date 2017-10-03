function v = VectorOfStressStrainMatrixInVM(m, varargin)
% VectorOfStressStrainMatrixInVM - Converts stress (strain) matrix to a
% vector in modifed Voigt-Mandel notation
%
%   USAGE:
%
%   t = VectorOfStressStrainMatrixInVM(v)
%
%   INPUT:
%
%   m
%       stress (strain) vector.
%
%   OUTPUT:
%
%   v
%       modified Voight-Mandel stress (strain) vector in [11, 22, 33, 23,
%       13, 12] order. Shears have factor of sqrt(2) multiplied to them.
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

if (size(m,1) == 3) && (size(m,2) == 3)
    s2  = sqrt(2);
    if strcmp(opts.Order, '11-22-33-23-13-12')
        v   = [ ...
            m(1,1); ...
            m(2,2); ...
            m(3,3); ...
            s2*m(2,3); ...
            s2*m(1,3); ...
            s2*m(1,2); ...
            ];
    elseif strcmp(opts.Order, '11-22-33-12-13-23')
        v   = [ ...
            m(1,1); ...
            m(2,2); ...
            m(3,3); ...
            s2*m(1,2); ...
            s2*m(1,3); ...
            s2*m(2,3); ...
            ];
    else
        disp('order not supported');
        return
    end
else
    disp('input matrix size incorrect');
end