function pout = pkfit_MapResult(pkpars, p)
% pkfitResultMapping - peak fit result mapping
%
%   USAGE:
%
%   pout = pkfitResultMapping(pkpars, p)
%
%   INPUT:
% 
%   pkpars
%       structure array with the following fields
% 
%       xdata       : data in the x-axis
%       pfunc_type  : peak function type (options are splitpseudovoigt,
%           gaussian, lorentzian, and pseudovoigt)
%       pbkg_order  : back ground order
%
%   p
%       peak fit result based on peak profile function
%
%   OUTPUT:
%
%   pout
%       Peak fit result (p) mapped to standard shape so that it can be
%       stored. For Gaussian and Lorentzian, the
%       mixing parameters are set to nan.
%
%   NOTE:
%   
%   Only works for a single peak for now.

pfunc_type  = pkpars.pfunc_type;

n   = pfunc_NumberOfParameters(pfunc_type);
switch lower(pfunc_type)
    case 'splitpseudovoigt'
        pout    = p;
    case 'gaussian'
        pout(1) = p(1);
        pout(2) = p(2);
        pout(3) = p(2);
        pout(4) = nan;
        pout(5) = nan;
        pout(6) = p(3);
        pout    = [pout; p((n+1):end)];
    case 'lorentzian'
        pout(1) = p(1);
        pout(2) = p(2);
        pout(3) = p(2);
        pout(4) = nan;
        pout(5) = nan;
        pout(6) = p(3);
        pout    = [pout; p((n+1):end)];
    case 'pseudovoigt'
        pout(1) = p(1);
        pout(2) = p(2);
        pout(3) = p(2);
        pout(4) = p(3);
        pout(5) = p(3);
        pout(6) = p(4);
        pout    = [pout; p((n+1):end)];
    otherwise
        disp('Unknown peak function!!')
end
