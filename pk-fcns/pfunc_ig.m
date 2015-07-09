function [p0, pLB, pUB] = pfunc_ig(params, ydata)
% pfunc_ig - generates initial guess for the peak fitting
%
%   USAGE:
%
%   f   = pkGaussian(p, params)
%
%   INPUT:
%
%   params
%       structure array with the following fields
% 
%       xdata       : data in the x-axis
%       pfunc_type  : peak function type (options are splitpseudovoigt,
%           gaussian, lorentzian, and pseudovoigt)
%       pbkg_order  : back ground order
%
%   ydata
%       peak intensity data corresponding to params.xdata
%
%   OUTPUT:
%
%   p0
%       initial guess
%
%   pLB
%       lower bound for the peak parameters
% 
%   pUB
%       upper bound for the peak parameters
%

pfunc_type  = params.pfunc_type;
pbkg_order  = params.pbkg_order;
[A, idx]    = max(ydata);

n       = pfunc_NumberOfParameters(pfunc_type);
xdata	= params.xdata;
xdata   = xdata(:);

pk_position = xdata(idx);

switch lower(pfunc_type)
    case 'splitpseudovoigt'
        numpk   = (length(p) - pbkg_order)/n;
        for i = 1:1:numpk
            ji  = n*(i - 1) + 1;
            jf  = n*i;
            
            ppk = p(ji:jf);
            ypk = pksplitpseudoVoigt(ppk,xdata);
            y   = y + ypk;
        end
        pbkg    = p((numpk*n+1):end);
        ybkg    = polyval(pbkg,xdata);
        
        y   = y + ybkg;
    case 'gaussian'
        numpk   = (length(p) - pbkg_order)/n;
        for i = 1:1:numpk
            ji  = n*(i - 1) + 1;
            jf  = n*i;
            
            ppk = p(ji:jf);
            ypk = pkGaussian(ppk,xdata);
            y   = y + ypk;
        end
        pbkg    = p((numpk*n+1):end);
        ybkg    = polyval(pbkg,xdata);
        
        y   = y + ybkg;
    case 'lorentzian'
        numpk   = (length(p) - pbkg_order)/n;
        for i = 1:1:numpk
            ji  = n*(i - 1) + 1;
            jf  = n*i;
            
            ppk = p(ji:jf);
            ypk = pkLorentzian(ppk,xdata);
            y   = y + ypk;
        end
        pbkg    = p((numpk*n+1):end);
        ybkg    = polyval(pbkg,xdata);
        
        y   = y + ybkg;
    case 'pseudovoigt'
        numpk   = (length(p) - pbkg_order)/n;
        for i = 1:1:numpk
            ji  = n*(i - 1) + 1;
            jf  = n*i;
            
            ppk = p(ji:jf);
            ypk = pkpseudoVoigt(ppk,xdata);
            y   = y + ypk;
        end
        pbkg    = p((numpk*n+1):end);
        ybkg    = polyval(pbkg,xdata);
        
        y   = y + ybkg;
    otherwise
        disp('Unknown peak function!!')
end