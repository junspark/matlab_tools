function fsym = SymmHarmonic(pb, cf, pts, sym)
% SYMMHARMONIC - Compute symmetrized harmonic.
%   
%   VERSION:  $Id: SymmHarmonic.m 170 2010-03-01 00:29:25Z boyce $
%
%   STATUS:  in development
%
%   USAGE:
%
%   fsym = SymmHarmonic(pb, cf, pts, sym)
%
%   INPUT:
%
%   pb is d x n
%      the result of PolyBasis()
%   cf is the coefficients
%   pts is a list of points to evaluate
%   sym is 4 x nsym, the crystal symmetries
%
%   OUTPUT:
%
%   fsym is 1 x npts
%        the symmetrized function
%
%   NOTES:
%
npts = size(pts, 2);
nfun = size(cf, 2);
%
fsym = zeros(npts, nfun);
printEvery = 10;
for i=1:size(cf, 2)
  if mod(i, printEvery) == 0
    fprintf('i = %d\n', i);
    if  (i ./ printEvery) == 10
      printEvery = 10 * printEvery;
    end
  end
  fsym(:, i) = SymmetrizeOne(pb, cf(:, i), pts, sym);
end
%
% ======================================== Helper functions
%
function h = SymmetrizeOne(pb, cf, pts, sym)
% SYMMETRIZEONE - 
%   
nsym = size(sym, 2);
npts = size(pts, 2);
%
hall = zeros(nsym, npts);
for i=1:nsym
  qpts = QuatProd(pts, repmat(sym(:, i), [1  npts]));
  hall(i, :) = PolyBasisEval(pb, cf, qpts);
end
h = mean(hall);
h(abs(h) < eps) = 0;
