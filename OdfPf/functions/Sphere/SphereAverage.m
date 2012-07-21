function [avg, optdat] = SphereAverage(pts, Pzation, x0, nlopts)
% SphereAverage - Find `average' of list of points on sphere.
%   
%   USAGE:
%
%   avg = SphereAverage(pts)
%   [avg, optdat] = SphereAverage(pts, Pzation, nlopts)
%
%   INPUT:
%
%   pts     is m x n, 
%           a list of n points in R^m of unit length
%   Pzation is a function handle,
%           (see note below)
%   x0      is the initial guess,
%           in given parameterization
%   nlopts  are options to be passed to the nonlinear minimizer.
%
%   OUTPUT:
%
%   avg is m x 1, 
%          is a unit vector representing the "average" of `pts'
%   optdat is a cell array, 
%          with three members,  {fval, exitflag, output}
%          (see documentation for `fminunc')
%
%   NOTES:
%
%   *  If only one argument is given, the average returned is 
%      the arithmetic average of the points.  If all three arguments
%      are given, then the average is computed using unconstrained
%      minimization of the sum of squared angles from the data points,
%      using the parameterization specified and the options given.
%      
%   *  See the matlab builtin `fminunc' for details.
%
%   *  This routine needs to be fixed.  Currently it uses the
%      parameterization given by `SpherePZ' instead of the 
%      function handle `PZation'.
%
avg = UnitVector(sum(pts, 2)); % arithmetic average
if (nargin == 1)
  return
end
%
%  Use nonlinear minimization.
%
fun = @SphDistanceFunc;
%
[x, fval, exitflag, output]  = fminunc(fun,x0, nlopts, pts, @SpherePZ);
optdat = {fval, exitflag, output};
%
avg = SpherePZ(x);
