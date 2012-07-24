function Re = ErrorRe(yData, yCalc, varargin)
% ErrorRe - calculates residual of the pattern.  for more details
% look at Rietvald Method (Young).
%
%   USAGE:
%
%   Re  = ErrorRe(yData, yCalc)
%
%   INPUT:
%   yData
%       experimental intensity data
%
%   yCalc
%       calculated intensity
%
%   OUTPUT:
%
%   Re
%       error

% default options
optcell = {...
    'Threshold', 0.0001, ...
    };

% update option
opts    = OptArgs(optcell, varargin);

i1  = yData > opts.Threshold;
i2  = yData < -opts.Threshold ;
i   = i1 | i2;

if sum(i) == 0
    Re  = 0;
else
    Re  = abs(yData(i) - yCalc(i))./abs(yData(i));
    Re  = mean(Re);
end