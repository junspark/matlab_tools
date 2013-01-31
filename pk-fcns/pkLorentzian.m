function [f, Jac] = pkLorentzian(p, x)
% pkLorentzian - symmetric peak profile function Lorentzian.  Refer to Rietveld
% Analysis text (Young) for more details.  The function is of the following
% form.
%
%   f = A*(c1^.5/Gamma/pi)./(1 + c1*x^2);
%   c0    : constant 4
%   A     : intensity
%   x     : (tth-tth_peak)/Gamma
%   Gamma : FWHM
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

c1      = 4;
A       = p(1);
Gamma   = p(2);
xPeak   = p(3);

delx    = (x - xPeak)./Gamma;

f       = (A*c1^.5/Gamma/pi)./(1+c1*delx.*delx);

% compute jacobian
% take the partial derivatives
dfdA        = f./A;

dfdGamma    = -A*c1^.5/Gamma/Gamma/pi./(1 + (c1*(x-xPeak).*(x-xPeak)./Gamma./Gamma)) + ...
    2*A*c1^(3/2)/Gamma/Gamma/Gamma/Gamma/pi.*(x-xPeak).*(x-xPeak)./(1 + (c1*(x-xPeak).*(x-xPeak)./Gamma./Gamma))./(1 + (c1*(x-xPeak).*(x-xPeak)./Gamma./Gamma));
% keyboard

dfdxPeak    = 2*A*c1^(3/2)/Gamma/Gamma/Gamma/pi.*(x-xPeak)./(1 + (c1*(x-xPeak).*(x-xPeak)./Gamma./Gamma))./(1 + (c1*(x-xPeak).*(x-xPeak)./Gamma./Gamma));

% arrange in Jacobian matrix
Jac = [...
    dfdA, ...
    dfdGamma, ...
    dfdxPeak...
    ];