function [f, Jac] = pkpseudoVoigt(p, x)
% pkpseudoVoigt - peak profile function Lorentzian.  Refer to Rietveld
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

A       = p(1);
Gamma   = p(2);
eta     = p(3);
xPeak   = p(4);

pG  = [A Gamma xPeak];
pL  = [A Gamma xPeak];

[fG, JacG]  = pkGaussian(pG, x);
[fL, JacL]  = pkLorentzian(pL, x);

f   = eta*fG + (1-eta)*fL;

% compute jacobian
dfdA        = eta*JacG(:,1) + (1-eta)*JacL(:,1);
dfdGamma    = eta*JacG(:,2) + (1-eta)*JacL(:,2);
dfdeta      = fG - fL;
dfdxPeak    = eta*JacG(:,3) + (1-eta)*JacL(:,3);

% arrange in Jacobian matrix
Jac = [...
    dfdA, ...
    dfdGamma, ...
    dfdeta, ...
    dfdxPeak, ...
    ];