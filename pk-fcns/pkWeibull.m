function f = pkWeibull(p, x)
% pkWeibull - Weibull distribution.  Refer to the Wikipedia page for more
% details. The function is of the following form.
%
%   f = (k/lambda)*((x - xpk)/lambda)^(k-1)e^(-(x-xpk)/lambda)^k) (x > xpk)
%
%   USAGE:
%
%   f = pkWeibull(p, x)
%
%   INPUT:
%   p
%       parameters for a peak.
%       p(1): shape parameter (k)
%       p(2): scale parameter (lambda)
%       p(3): peak position (xpk)
%
%   x
%       x coordinates typically in 2 theta, mm, or pixels
%
%   OUTPUT:
%
%   f
%       value of peak function at each x
%

p   = p(:);
x   = x(:);

k       = p(1);
lambda  = p(2);
xpk     = p(3);

delx    = (x - xpk)./lambda;

f       = (k / lambda)*((delx).^(k - 1)).* exp( -delx.^k);