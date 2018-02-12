function [F, J] = gs1(p,x)
% [F, J] = gs1(p,x):  Gauss distribution
% INPUT
% p: FLOAT(3), parameter values
%	p(1): peak value
%	p(2): center position
%	p(3): FWHM
% x: FLOAT(npts,1), x data points
% OUTPUT
% F: FLOAT(npts), function values
% J: FLOAT(npts,3), partial derivatives
% COMMENTS
% Can be called from SUMFUN1 for LSQ fitting

x = x(:);
F = p(1)*exp(-4*log(2)*((x-p(2))/p(3)).^2);

if nargout == 2 % Calculate derivatives only when requested
   arg1 = ((x-p(2))/p(3)).^2;
   arg2 = exp(-4*log(2)*arg1);
   J(:,1) = arg2; 
   J(:,2) = ( 8*log(2)*p(1)*(x-p(2))/p(3)^2.*arg2 );
   J(:,3) = ( 8*log(2)*p(1)*arg1/p(3).*arg2 );
end