function [func, Jac] = pfunc_Lorentzian(p, x, jacFlag)
% symmetric peak profile function Gaussian
% refer to Rietveld book for more details
% p is the parameters for a peak
% x is the x axis (typically in radial distance or in tth)
% jacFlag is the Jacobian flag
%
% Lorentzian peak function
% f = A*(c1^.5/Gamma/pi)./(1 + c1*x^2);
%
% c0    : constant 4*log(2)
% A     : intensity
% x     : (tth-tth_peak)/Gamma
% Gamma : FWHM

p   = p(:);
x   = x(:);

c1      = 4;

A       = p(1);
Gamma   = p(2);
xPeak   = p(3);

delx    = (x - xPeak)./Gamma;

func    = A*c1^.5/Gamma/pi./(1+c1*delx.*delx);

% compute jacobian
if nargout > 1

    % take the partial derivatives
    dfdA        = func./A;
    dfdGamma    = -func./Gamma + ...
        2*c1^.5*pi/A*delx.*delx.*func.*func;
    dfdxPeak    = (2*c1^.5*delx*pi/A).*func.*func;

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
