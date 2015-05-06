function pout = pkfitResultMapping(pkpars, p)
% pkfitResultMapping - peak fit result mapping
%
%   USAGE:
%
%   pout = pkfitResultMapping(pkpars, p)
%
%   INPUT:
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
%       peak fit result (p) mapped to standard shape so that it can be
%       stored. for splitpseudoVoigt, the two widths and mixing parameters
%       are averaged. for gaussian and lorentzian, mixing parameters are
%       set to zeros.
%
%   NOTE:
%   
%   Only works for a single peak for now.

pfunc_type  = pkpars.pfunc_type;

switch lower(pfunc_type)
    case 'splitpseudovoigt'
        pout(1) = pr(1);
        pout(2) = mean(p(2:3));
        pout(3) = mean(p(4:5));
        pout(4) = pr(6);
        pout    = [pout; pr(7:end)];
    case 'gaussian'
        pout(1) = pr(1);
        pout(2) = 0;
        pout(3) = pr(2);
        pout(4) = pr(3);
        pout    = [pout; pr(4:end)];
    case 'lorentzian'
        pout(1) = pr(1);
        pout(2) = 0;
        pout(3) = pr(2);
        pout(4) = pr(3);
        pout    = [pout; pr(4:end)];
    case 'pseudovoigt'
        pout    = p;
    otherwise
        disp('Unknown peak function!!')
end
