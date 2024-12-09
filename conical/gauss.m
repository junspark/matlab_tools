function y = gauss(p, x)

%Gaussian
y = p(1) .*exp(-log(2)*4*((x-p(2)) ./p(3)) .^2) + p(4) + p(5).*(x-p(2));