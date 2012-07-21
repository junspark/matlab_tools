function all = ToAllNodes(red, eqv)
% ToAllNodes - Spread array at independent nodes to all nodes.
%   
%   USAGE:
%
%   all = ToAllNodes(red, eqv)
%
%   INPUT:
%
%   red is m x n, 
%       a list of n m-vectors at the independent nodes of a mesh
%
%   eqv is 2 x l, (integer)
%       the list of node equivalences (new #, old #)
%
%   OUTPUT:
%
%   all is m x k, 
%       the array at all nodes of the mesh; k = n + l
%
if (isempty(eqv))
  all = red;
  return
end
%
%  If the input is simply a vector, then preserve shape on output.
%
[m n] = size(red);
if (n == 1) 
  transpose = 1;
  red = red';
  [m n] = size(red);
else
  transpose = 0;
end
%
eqvmax = max(eqv(1, :)); 
eqvmin = min(eqv(1, :)); 
%
if ( (n+1) ~= eqvmin)
  error('Data array does not match equivalence array')
end
%
nall = eqvmax;
%
all               = zeros(m, nall);    % allocate odf
all(:, 1:n)       = red;
all(:, eqv(1, :)) = red(eqv(2, :));
%
if (transpose)       % return vector of same type (row/col)
  all = all'; 
end
