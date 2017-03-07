function y = pkLognormal(p, x)
% pkLognormal - Log normal function. 
%
%   USAGE:
%
%   y   = pkLognormal(p, x)
%
%   INPUT:
%   p
%       parameters for a peak.
%       p(1): amplitude
%       p(2): log mean
%       p(3): log standard deviation
%
%   x
%       x coordinates typically in 2theta, mm, or pixels
%
%   OUTPUT:
%
%   f
%       value of peak function at each x
%

A   = p(1);
mu  = p(2);
sig = p(3);

idx = x < 0;
idx = sum(idx);
if idx > 0
    disp('x > 0');
    return
end

a   = A./x;
a   = a./sig;
a   = a./sqrt(2*pi);

b   = log(x);
b   = b - mu;
b   = -b.*b;
b   = b./2;
b   = b./(sig*sig);

y   = a.*exp(b);
