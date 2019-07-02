function q = d2q(d)
% d2q - converts d-spacings to q
%
%   USAGE:
%
%   q = d2Q(d)
%
%   INPUT:
%
%   d
%       d-spacing (length)
%
%   OUTPUT:
%
%   q 
%       momentum transfer (length^-1)

d   = d(:);
q   = (2*pi)./d;