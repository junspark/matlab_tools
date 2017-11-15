function [F, J, area] = psv1(p,x)
% [F, J] = psv1(p,x):  PseudoVoigt distribution
% INPUT
% p: DOUBLE(4), parameter values
%	p(1): peak value
%	p(2): center position
%	p(3): FWHM
%  p(4): profile shape factor (0-1:Gaussian-Lorentzian)
% x: DOUBLE(npts,1), x data points
% OUTPUT
% F: DOUBLE(npts), function values
% J: Double(npts,4), partial derivatives
% area: DOUBLE, peak area
% COMMENTS
% Can be called from SUMFUN1 for LSQ fitting
%
% UL 20010118

x = x(:);
g = exp(-.5*((x-p(2))*2*sqrt(2*log(2))/p(3)).^2);
l = 1 ./ (1 + ((x - p(2))*2/p(3)).^2);
F = p(1) * (p(4)*l + (1 - p(4))*g);

if nargout >= 1 % Calculate derivatives only when requested
  c1	= 2*sqrt(2*log(2));
  arg1	= 2*(x-p(2))/p(3);
  arg2	= c1*(x-p(2))/p(3);
  l	 = 1 + arg1.^2;
  g	 = exp(-.5*arg2.^2);
  J(:,1) = p(4)./l + (1-p(4))*g;
  J(:,2) = p(1)*(p(4)*l.^(-2)*4.*arg1/p(3) + (1-p(4))*g.*arg2*c1/p(3));
  J(:,3) = (x-p(2))/p(3).*J(:,2);
  J(:,4) = p(1)*(1./l - g);
  end

if nargout >= 2
  area	= p(1)*p(3)*((1 - p(4))*sqrt(pi/log(2))/2 + p(4)*pi/2);
  end
