function hsym = SymHKL(hkl, sym, invsym)
% SYMHKL - Find symmetric equivalents of the vector hkl.
%
%   hsym = SymHKL(hkl, sym)
%   hsym = SymHKL(hkl, sym, invsym)
%
%   hkl is 3 x 1
%	a vector
%   sym is 4 x n
%	the symmetry group in quaternions
%   invsym is a scalar (optional)
%       if present and nonzero, it includes inversion 
%       symmetries (h --> -h)
%
%   hsym is 3 x m
%	the list of symmetrically equivalent vectors to hkl,
%       with duplicates removed;  vectors are compared using
%	the UniqueVectors function to a tolerance of 1.0e-14
%   
nsym   = size(sym, 2);
Rsym   = RMatOfQuat(sym);

allhkl = MultMatArray(Rsym, repmat(hkl, [1 1 nsym]));
allhkl = squeeze(allhkl);
%
if nargin < 3
  invsym = 0;
end
%
if invsym == 0
  hsym   = UniqueVectors([allhkl]);
else
  hsym   = UniqueVectors([allhkl, -allhkl]);
end
