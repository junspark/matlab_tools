function [cfs, pb] = PolyHarmonicBasis(dim, deg)
% PolyHarmonicBasis - Basis for spherical harmonics
%   
%   VERSION:  $Id: PolyHarmonicBasis.m 170 2010-03-01 00:29:25Z boyce $
%
%   STATUS:  in development
%
%   USAGE:
%
%   cfs = PolyHarmonicBasis(dim, deg)
%   [cfs, pb] = PolyHarmonicBasis(dim, deg)
%
%   INPUT:
%
%   dim is a scalar integer,
%	the dimension of the underlying Euclidean space
%   deg is a scalar integer,
%       the degree of the polynomial space
%
%   OUTPUT:
%
%   cfs is n x m
%       It is a set of coefficients of the PolyBasis functions
%       which produce zero LaPlacian.
%
%   pb  is dim x n (optional)
%       the polynomial basis of dimension dim and degree deg.
%
%   NOTES:
%
pb  = PolyBasis(dim, deg);
lp  = PolyLaPlacian(pb);
cfs = null(full(lp));
