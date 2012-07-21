function der = RefDerivatives(n)
% RefDerivatives - Reference shape function derivatives.
%   
%   USAGE:
%
%   der = RefDerivatives(n)
%
%   INPUT:
%
%   n is a positive integer,
%     the dimension of the simplex
%
%   OUTPUT:
%
%   der is n x (n+1),
%	column j contains the gradient of the j'th
%       barycentric coordinate with respect to the 
%       coordinate directions
%
der = [eye(n), repmat([-1], [n 1])];
