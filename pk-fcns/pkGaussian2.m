function f = pkGaussian2(p, r_eta)
% pkGaussian2 - symmetric peak profile function Gaussian in 2D.  
%
%   f = A * exp( -(r-r0)^2 / 2/gamma_r^2 -(eta-eta0)^2 / 2/gamma_eta^2) + B
%
%   USAGE:
%
%   f = pkGaussian(p, x)
%
%   INPUT:
%   p
%       parameters for a peak.
%       p(1): amplitude
%       p(2): gamma_r (~width, FWHM in radial direction)
%       p(3): gamma_eta (~width, FWHM in azimuthal direction)
%       p(4): r0 (peak position in radial direction)
%       p(5): eta0 (peak position in eta direction)
%       p(6): B (constant background)
%
%   r_eta
%       r and eta coordinates [n x 2 array: column 1 is r and column 2 is eta]
%
%   OUTPUT:
%
%   f
%       value of peak function at each x

p   = p(:);
r   = r_eta(:,1);
eta = r_eta(:,2);

A           = p(1);
Gamma_r     = p(2);
Gamma_eta   = p(3);
rPeak       = p(4);
etaPeak     = p(5);
B           = p(6);

delr    = (r - rPeak)./Gamma_r;
deleta  = (eta - etaPeak)./Gamma_eta;

f       = A*exp(-0.5*delr.*delr + -0.5*deleta.*deleta) + B;