function [dh, eVals] = DiscreteHarmonics(m, nHarm)
% DISCRETEHARMONICS - Generate the first several discrete harmonics.
%   
%   VERSION:  $Id: DiscreteHarmonics.m 170 2010-03-01 00:29:25Z boyce $
%
%   STATUS:  in development
%
%   USAGE:
%
%   [dh, eVals] = DiscreteHarmonics(m, nHarm)
%
%   INPUT:
%
%   m is a MeshStructure
%        on either the sphere S^2, the sphere S^3 or Rodrigues space;
%        the mesh is expected to have 'h1form' and 'l2ip' matrices as fields
%
%   nHarm is an integer
%            the number of requested harmonics
%
%   OUTPUT:
%
%   dh is n x nHarm
%         the discrete harmonics
%
%   NOTES:
%
%   *  Calling 'eigs' with opts.v0 set to ones (as here)  provides a consistent answer here;
%      Normally, v0 is randomly chosen, and so it is possible to get different
%      bases, especially given that the abstract problem has great multiplicity.
%
thisFile = mfilename();
%
%  Options:  
%
OPTS.disp   = 0;               %  display level:  0
OPTS.v0     = ones(m.numind, 1); %  initial guess for consistency

% enforce symmetric h1form and l2ip matrices
m.h1form    = 0.5*(m.h1form + m.h1form');
m.l2ip      = 0.5*(m.l2ip + m.l2ip');

% compute shm
disp('computing eigenvalues ...'); tic
[dh, D] = eigs(m.h1form, m.l2ip, nHarm, 'SM', OPTS);

[eVals, idx]    = sort(diag(D));
dh = dh(:, idx); % reverse order on output
disp('done'); toc