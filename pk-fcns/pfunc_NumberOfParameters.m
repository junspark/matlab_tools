function NumberOfParameters = pfunc_NumberOfParameters(pfunc_type)
% pfunc_NumberOfParameters - outputs the number of parameters used in a peak profile
% function
%
%   USAGE:
%
%   nPrms = pfunc_NumberOfParameters(PkFcn)
%
%   INPUT:
%
%   pfunc_type
%       peak profile function name (string). Options are splitpseudovoigt,
%       gaussian, lorentzian, lognormal, and pseudovoigt.
% 
%   OUTPUT:
%
%   NumberOfParameters
%       number of parameters
%

switch lower(pfunc_type)
    case 'splitpseudovoigt'
        NumberOfParameters   = 6;
    case 'pseudovoigt'
        NumberOfParameters   = 4;
    case 'gaussian'
        NumberOfParameters   = 3;
    case 'lorentzian'
        NumberOfParameters   = 3;
    case 'lognormal'
        NumberOfParameters   = 3;
    otherwise
        disp('Unknown peak function!!')
end