function [sk, gk, Hk] = SpherePZ(x)
% SpherePZ - Generate point, gradients and Hessians of map to sphere.
%   
%   USAGE:
%
%   sk           = SpherePZ(x)
%   [sk, gk]     = SpherePZ(x)
%   [sk, gk, Hk] = SpherePZ(x)
%
%   INPUT:
%
%   x is d x 1, 
%     a vector with norm <= 1
%
%   OUTPUT:
%
%   sk is e x 1, 
%      a point on the sphere (sqrt(1-x^2), x)
%   gk is d x e, 
%      the gradients of each component of sk
%   Hk is d x d x e, 
%      the Hessians of each component of sk
%
x = x(:);         % make column vector
d = length(x);
%
%  Components.
%
Nx = sqrt(1 - x'*x);
sk = [Nx; x];
%
if (nargout == 1)
  return
end
%
%  Gradients of components.
%
mNxi = -1/Nx;
gN   = mNxi * x; 
gk   = [gN eye(d)];
%
if (nargout == 2)
  return
end
%
%  Hessians of components.
%
Hk = [mNxi*(eye(d) + gN*gN') zeros(d, d*d)];
Hk = reshape(Hk, [d d (d+1)]);
%
return
