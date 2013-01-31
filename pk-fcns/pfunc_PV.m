function func = pfunc_PV(p, x)
% symmetric peak profile function pseudo Voigt
% refer to Rietveld book for more details
% p is the parameters for a peak
% x is the x axis (typically in radial distance or in tth)
%
% pseudo-Voigt peak function
% f = A*(1-n)*G + A*n*L
% where G and L are Gaussian and Lorentzian functions sharing the same peak
% position and amplitude
%
% c0    : constant 4*log(2)
% A     : intensity
% x     : (tth-tth_peak)/Gamma
% Gamma : FWHM

p = p(:);
x = x(:);

A  = p(1);
n  = p(2);
G  = p(3);
x0 = p(4);

% the TCH pseudo-Voight:
%
% f = A*((1 - n)*fg + n*fl)

sig2fwhm = 2*sqrt(log(2));
gam2fwhm = 2;

gamg = G/sig2fwhm;
gaml = G/gam2fwhm;

delx = x - x0;

Ag = 1/sqrt(pi);
Al = 1/pi;

fg = Ag/gamg*exp(-delx.^2/gamg^2);
fl = Al/gaml ./ (1 + (delx/gaml).^2);

func = A*( (1 - n)*fg + n*fl );
