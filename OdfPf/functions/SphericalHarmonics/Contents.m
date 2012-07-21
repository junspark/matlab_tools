%
%  SPHERICAL HARMONICS.
%
%  This directory contains functions for evaluating spherical
%  harmonics either analytically in terms of hnomgeneous polynomials 
%  or discretely through eigenvectors on meshes over orientation space.
%
%  ANALYTIC EVALUATION ROUTINES.
%
% The main routines in this section is PolyHarmonicBasis and SymmHarmonic.
%
% PolyBasis(dim, deg) -- create a list multi-indices (exponents) spanning polynomials
%          -- of given dimension and degree
% PolyBasisDim(dim, deg) -- return dimension of space of homogeneous polynomials
% PolyBasisIndex(mi) -- return the linear (scalar) index of the input multi-index
% PolyLaPlacian(mi) -- computes LaPlacian for an entire space
% PolyBasisEval(pbasis, coeffs, pts) -- evaluate the polynomial at given points, where
%          -- the polynomial is determined by the given coefficients on the given basis
%
% PolyHarmonicBasis(dim, deg) -- returns coefficients and optionally polynomial basis 
%          -- for harmonic polynomials
% SymmHarmonic(pb, cf, pts, sym) -- symmetrizes harmonics 
%
%  DISCRETE HARMONIC ROUTINES
% 
% DiscreteHarmonics()
%
% ExportDX1()
% ExportDXField()
% MakeHarmonics()
% SphL2IP()
% 
% WriteDXHarmonics()
