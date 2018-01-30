function [F, J] = backg1(p,x)

x = x(:) - (x(1)+ x(length(x)))/2;
F = zeros(size(x));
np = length(p);
for i=1:np
	F = F + p(np+1-i) * x.^(i-1);
	if nargout == 2 %  F and J requested
	J(:,i) = x.^(np-i);
	end
end

%F=F+abs(p(1)*length(x)/2); %JA added Feb05 to better match background  % removed by AC Mar11

%if nargin == 4 % put out weighted difference
%	y = y(:); w = w(:);
%	F = (F - y)./w;
%	J = J ./ (w*ones(1,np));
%end
