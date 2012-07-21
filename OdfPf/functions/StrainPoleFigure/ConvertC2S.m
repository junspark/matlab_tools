function s = ConvertC2S(c, varargin)
% ConvertC2S - Converts elastic stiffness constants to elastic compliance
% constants
%   
%   USAGE:
%
%   s   = ConvertC2S(c)
%   s   = ConvertC2S(c, 'Symmetry', symmetry group)
%
%   INPUT:
%
%   c   = (3 x 1) vector with [c11 c12 c44];
%
%   OUTPUT:
%
%   s   = (3 x 1) vector with [s11 s12 s44];
%


% default options
optcell = {...
    'Symmetry', 'Cubic', ...
    };

% update option
opts    = OptArgs(optcell, varargin);

if strcmp(lower(opts.Symmetry), 'cubic')
    c11 = c(1);
    c12 = c(2);
    c44 = c(3);
    
    d   = (c11*c11 + c11*c12 - 2*c12*c12);
    
    s11 = (c11 + c12) / d;
    s12 = -c12 / d;
    s44 = 1/c44;
    
    s   = [s11; s12; s44];
else
    disp('symmetry group not supported')
end