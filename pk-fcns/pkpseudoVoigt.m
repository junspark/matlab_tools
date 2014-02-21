function [f, Jac] = pkpseudoVoigt(p, x)
% pkpseudoVoigt - peak profile function pseudo-Voigt.  Refer to Rietveld
% Analysis text (Young) for more details. This is a mix of Guassian and
% Lorentzian functions.
%
%   USAGE:
%
%   [f, Jac] = pkpseudoVoigt(p, x)
%
%   INPUT:
%   p
%       parameters for a peak.
%       p(1): amplitude
%       p(2): gamma (~width, FWHM)
%       p(3): mixing parameter between Gaussian and Lorentzian
%       P(4): peak position
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