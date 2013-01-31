function [f, Jac] = pkGaussian(p, x)
% pkGaussian - symmetric peak profile function Gaussian.  Refer to Rietveld
% Analysis text (Young) for more details.  The function is of the following
% form.
%
%   f = A*(c0^.5/Gamma/pi^.5).*exp(-c0*x.*x)
%
%   USAGE:
%
%   [f, Jac] = pkGaussian(p, x)
%
%   INPUT:
%   p
%       parameters for a peak.
%       p(1): amplitude
%       p(2): gamma (~width, FWHM)
%       P(3): peak position
%
%   x
%       x coordinates typically in 2 theta, mm, or pixels
%
%   OUTPUT:
%
%   f
%       value of peak function at each x
%
%   Jac
%       Jacobian of the peak function

p   = p(:);
x   = x(:);

c0      = 4*log(2);
A       = p(1);
Gamma   = p(2);
xPeak   = p(3);

delx    = (x - xPeak)./Gamma;

f       = (A/Gamma)*(c0/pi)^.5*exp(-c0*delx.*delx);

% compute jacobian
% take the partial derivatives
dfdA        = f./A;
dfdGamma    = -f./Gamma + ...
    (2*c0/Gamma*delx.^2).*f;
dfdxPeak    = (2*c0/Gamma*delx).*f;

% arrange in Jacobian matrix
Jac = [...
    dfdA, ...
    dfdGamma, ...
    dfdxPeak...
    ];