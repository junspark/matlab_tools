function f = pksplitpseudoVoigt(p, x)
% pksplitpseudoVoigt - peak profile function split pseudo-Voigt. This is a
% mix of Gaussian and Lorentzian functions. 
%
%   USAGE:
%
%   f   = pkpseudoVoigt(p, x)
%
%   INPUT:
%   p
%       parameters for a peak.
%       p(1): amplitude
%       p(2): gamma left hand side (~width, FWHM)
%       p(3): gamma right hand side (~width, FWHM)
%       p(4): mixing parameter between Gaussian and Lorentzian left hand side
%       p(5): mixing parameter between Gaussian and Lorentzian right hand side
%       P(6): peak position (demarcation for left and right and sides)
%
%   x
%       x coordinates typically in 2 theta, mm, or pixels
%
%   OUTPUT:
%
%   f
%       value of peak function at each x

p   = p(:);
x   = x(:);
f   = 0.*x;

c1          = 4;

A           = p(1);
Gamma_LHS   = p(2);
Gamma_RHS   = p(3);
eta_LHS     = p(4);
eta_RHS     = p(5);
xPeak       = p(6);

AG_LHS  = A*Gamma_LHS;
AG_RHS  = A*Gamma_RHS;

AL_LHS  = (A/c1^.5)*Gamma_LHS*pi;
AL_RHS  = (A/c1^.5)*Gamma_RHS*pi;

idx = find(x < xPeak);
idx = idx(end);

x_LHS   = x(1:idx);
x_RHS   = x((idx+1):end);

pG_LHS  = [AG_LHS Gamma_LHS xPeak];
pL_LHS  = [AL_LHS Gamma_LHS xPeak];
pG_RHS  = [AG_RHS Gamma_RHS xPeak];
pL_RHS  = [AL_RHS Gamma_RHS xPeak];

fG_LHS  = pkGaussian(pG_LHS, x_LHS);
fL_LHS	= pkLorentzian(pL_LHS, x_LHS);
fG_RHS	= pkGaussian(pG_RHS, x_RHS);
fL_RHS	= pkLorentzian(pL_RHS, x_RHS);

f_LHS   = eta_LHS*fG_LHS + (1-eta_LHS)*fL_LHS;
f_RHS   = eta_RHS*fG_RHS + (1-eta_RHS)*fL_RHS;

f(1:idx)        = f_LHS;
f((idx+1):end)  = f_RHS;