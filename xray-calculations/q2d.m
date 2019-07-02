function d = q2d(q)
% d2q - converts q to d-spacings
%
%   USAGE:
%
%   d = q2d(q)
%
%   INPUT:
%
%   q 
%       momentum transfer (length^-1)
%
%   OUTPUT:
%
%   d
%       d-spacing (length)

q   = q(:);
d   = (2*pi)./q;