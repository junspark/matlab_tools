function y = pfunc_switch(p, params)
% pfunc_switch - switch objective function for peak fitting
%
%   USAGE:
%
%   f   = pkGaussian(p, params)
%
%   INPUT:
%   p
%       parameters for a peak appropriate for the peak function of choice
%   
%   xr
%       data in the x-axis
%   
%   params
%       structure array with the following fields
% 
%       pfunc_type  : peak function type (options are splitpseudovoigt,
%           gaussian, lorentzian, and pseudovoigt)
%       pbkg_order  : back ground order
%
%   OUTPUT:
%
%   f
%       value of peak function at each x
%
p	= p(:);

xdata   = params.xdata;
% xdata   = xr;
xdata   = xdata(:);

y       = xdata.*0;

pfunc_type  = params.pfunc_type;
pbkg_order  = params.pbkg_order;

n   = pfunc_NumberOfParameters(pfunc_type);
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
    case 'lognormal'
        numpk   = (length(p) - pbkg_order)/n;
        for i = 1:1:numpk
            ji  = n*(i - 1) + 1;
            jf  = n*i;
            
            ppk = p(ji:jf);
            ypk = pkLognormal(ppk,xdata);
            y   = y + ypk;
        end
        pbkg    = p((numpk*n+1):end);
        ybkg    = polyval(pbkg,xdata);
        
        y   = y + ybkg;
    otherwise
        disp('Unknown peak function!!')
end