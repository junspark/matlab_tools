function f = pfunc_switch(p, params)
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
%   params
%       structure array with the following fields
% 
%       xdata       : data in the x-axis
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
xdata   = xdata(:);
f       = xdata.*0;

pfunc_type  = params.pfunc_type;
pbkg_order  = params.pbkg_order;

switch lower(pfunc_type)
    case 'splitpseudovoigt'
        numpk   = (length(p) - pbkg_order)/6;
        for i = 1:1:numpk
            ji  = 6*(i - 1) + 1;
            jf  = 6*i;
            
            ppk = p(ji:jf);
            ypk = pksplitpseudoVoigt(ppk,xdata);
            f   = f + ypk;
        end
        pbkg    = p((numpk*6+1):end);
        ybkg    = polyval(pbkg,xdata);
        
        f   = f + ybkg;
    case 'gaussian'
        numpk   = (length(p) - pbkg_order)/3;
        for i = 1:1:numpk
            ji  = 3*(i - 1) + 1;
            jf  = 3*i;
            
            ppk = p(ji:jf);
            ypk = pkGaussian(ppk,xdata);
            f   = f + ypk;
        end
        pbkg    = p((numpk*pbkg_order+1):end);
        ybkg    = polyval(pbkg,xdata);
        
        f   = f + ybkg;
    case 'lorentzian'
        numpk   = (length(p) - pbkg_order)/3;
        for i = 1:1:numpk
            ji  = 3*(i - 1) + 1;
            jf  = 3*i;
            
            ppk = p(ji:jf);
            ypk = pkLorentzian(ppk,xdata);
            f   = f + ypk;
        end
        pbkg    = p((numpk*pbkg_order+1):end);
        ybkg    = polyval(pbkg,xdata);
        
        f   = f + ybkg;
    case 'pseudovoigt'
        numpk   = (length(p) - pbkg_order)/4;
        for i = 1:1:numpk
            ji  = 4*(i - 1) + 1;
            jf  = 4*i;
            
            ppk = p(ji:jf);
            ypk = pkpseudoVoigt(ppk,xdata);
            f   = f + ypk;
        end
        pbkg    = p((numpk*pbkg_order+1):end);
        ybkg    = polyval(pbkg,xdata);
        
        f   = f + ybkg;
    otherwise
        disp('Unknown peak function!!')
end