function [func, Jac] = pfunc_Gaussian(p, x, jacFlag)
% symmetric peak profile function Gaussian
% refer to Rietveld book for more details
% p is the parameters for a peak
% x is the x axis (typically in radial distance or in tth)
% jacFlag is the Jacobian flag
%
% Gaussian peak function
% f = A*(c0^.5/Gamma/pi^.5).*exp(-c0*x.*x);
%
% c0    : constant 4*log(2)
% A     : intensity
% x     : (tth-tth_peak)/Gamma
% Gamma : FWHM

p   = p(:);
x   = x(:);

c0          = 4*log(2);
A           = p(1);
Gamma   = p(2);
xPeak   = p(3);

delx    = (x - xPeak)./Gamma;

func    = A*(c0^.5/Gamma/pi^.5)*exp(-c0*delx.*delx);

% compute jacobian
if nargout > 1

    % take the partial derivatives
    dfdA        = func./A;
    dfdGamma    = func./Gamma + ...
        (2*c0/Gamma*delx.^2).*func;
    dfdxPeak    = (2*c0/Gamma*delx).*func;

    % arrange in Jacobian matrix
    Jac = [...
        dfdA, ...
        dfdGamma, ...
        dfdxPeak...
        ];

    if nargin > 2
        if ~islogical(jacFlag)
            jacFlag = logical(jacFlag);
        end
        Jac = Jac(:, jacFlag);
    end
end
