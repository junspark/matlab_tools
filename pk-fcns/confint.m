function [conf,var]=confint(x,f,j,gamma)
%CONFINT Calculates confidence intervals of least squares estimation.
%
%	[CONF,VAR] = CONFINT(X,R,J) returns 95% confident interval CONF
%	and variance VAR for the least mean square estimation X, 
%	estimation residuals R and the Jacobian matrix J.
%
%	[CONF,VAR] = CONFINT(X,R,J,GAMMA) returns the confident interval and
%	variance with specified confident coefficient GAMMA (0 < GAMMA < 1).
%	The default value for GAMMA is 0.95.
%
%	The confident interval calculation is valid for "over determined"
%	estimation only. The estimated solution can be either linear or 
%	non-linear least square estimations.
%
%	Examples to estimate the confidence interval.
%	For a nonlinear least squares problem use:
%		[x,options,f,j]=leastsq('fun',x0,...);
%		[conf, var] = confint(x,f,j);
%
%	For a linear case:
%		A*x=b;
%		x=A\b;
%		[conf, var] = confint(x,b-A*x,A);
%
%	See also LEASTSQ.
%

%	Wes Wang, Andy Grace 6-28-93
%	Copyright (c) 1993 by the MathWorks, Inc.

%initialization
if nargin < 3
   disp('usage: [CONF,VAR] = CONFINT(X,S,J)');
elseif nargin <4
   gamma = .95;
end;
f = f(:);
[m,n] = size(j);
if m<n
   error('CONFINT is valid for over determined estimation problem only');
end;
if length(f) ~= m
   error('dimension R and J are not consistent')
end;
if length(x) ~= n
   error('dimension X and J are not consistent')
end;

% approximation when a column is zero vector
temp = find(max(abs(j))==0);
if ~isempty(temp)
   j(temp,:) = j(temp,:) + sqrt(eps);
end;

%calculate covariance
HESS = j' * j;
covar = f' * f * inv(HESS) / (m - n);
var = diag(covar);

% calculate student distribution
%reference statistic toolbox for the algorithm.
v = m-n;
temp = exp(gammaln((v + 1) / 2) - gammaln(v/2));
v = temp / (sqrt(v*pi) * (1 + (((1-gamma)/2)^ 2) / v)^ ((v + 1)/2));

% calculate confidence interval
temp = sqrt(var) .* v;
conf = [x(:) - temp x(:) + temp];

%--end of confint.m---

