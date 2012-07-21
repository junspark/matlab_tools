function [f, gf, Hf] = SphDistanceFunc(x, pts, Sofx)
% SphDistanceFunc - Return half sum of squared distances on sphere.
%   
%   USAGE:
%
%   f           = SphDistanceFunc(x, pts, @Sofx)
%   [f, gf]     = SphDistanceFunc(x, pts, @Sofx)
%   [f, gf, Hf] = SphDistanceFunc(x, pts, @Sofx)
%
%   INPUT:
%
%   x    is d x 1, 
%        a point in parameter space 
%   pts  is (d+1) x n, 
%        a list of n points on the sphere
%   Sofx is a function handle, 
%        returning parameterization component quantities 
%        (function, gradient and Hessian)
%
%   OUTPUT:
%
%   f  is a scalar, 
%      the objective function at x
%   gf is a vector, 
%      the gradient of f at x
%   Hf is a matrix, 
%      the Hessian of f at x
%
%   NOTES:
%
%   *  See MisorientationStats
%
global neval
neval = neval + 1;
%
[d1 n] = size(pts);
d      = d1 - 1;  % dimension of parameter space
%
%  Compute point according to parameterization.
%
if     (nargout == 1)          % function value only
  s           = feval(Sofx, x);
elseif (nargout == 2)          % and gradients
  [s, gs]     = feval(Sofx, x);
elseif (nargout == 3)          % and Hessians
  [s, gs, Hs] = feval(Sofx, x);
  Hs = reshape(Hs, [d*d (d+1)]);  % more efficient form for later use
end
%
%  Return function value.
%
ctheta = min(1, s'*pts);
thetai = acos(ctheta); % row vector
%
f = 0.5 * (thetai * thetai');          % function value
%
if (nargout == 1)
  return
end
%
%  Compute gradient.
%
gc = gs * pts;
stheta = sin(thetai);
%
limit  = (thetai <= eps);
%
thfac1(~limit) = thetai(~limit)./stheta(~limit);   
thfac1( limit) =  1;  % replace any NaN's generated above
%
gf = gc * thfac1';
gf = -gf;
%
if (nargout == 2)
  return
end
%
%  Compute Hessian.
%
Hc = Hs * pts;
%
limit  = (thetai <= eps^(1/3));
%
thfac3        = (stheta - (thetai.*ctheta))./ (stheta .^ 3);
thfac3(limit) = 1/3;   % replace NaN's
%
gcgct = reshape(RankOneMatrix(gc), [d*d n]);
%
Hf = gcgct * thfac3' - Hc * thfac1';
Hf = reshape(Hf, [d d]);
%
return
