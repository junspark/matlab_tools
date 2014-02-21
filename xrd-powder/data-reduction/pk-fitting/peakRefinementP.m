function funcOut = peakRefinementP(x,xdata)
% A=x(1);
% mix=x(2);
% G=x(3);
% x0=x(4);
% bkG1=x(5);
% bkG2=x(6);
bkGrndI = xdata.*x(5)+x(6);

%   p
%       parameters for a peak.
%       p(1): amplitude
%       p(2): gamma (~width, FWHM)
%       p(3): mixing parameter between Gaussian and Lorentzian
%       P(4): peak position

bkGrndI=xdata.*x(5)+x(6);

func    = pkpseudoVoigt([x(1) x(2) x(3) x(4)], xdata);
spectra = bkGrndI + func';

funcOut = spectra;
