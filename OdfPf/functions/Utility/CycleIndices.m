function cycle = CycleIndices(n)
% CycleIndices - Cycle the indices 1:n.
%   
%   USAGE:
%
%   cycle = CycleIndices(n)
%
%   INPUT:
%
%   n is a postive integer
%
%   OUTPUT:
%
%   cycle is an n x n array of integers ,
%         the columns are the indices 1:n cyclically permuted
%
ident = 1:n;
%
cols = repmat(ident, n, 1);
rows = cols';
%
cycle = 1 + mod(cols + rows - 2, n);
