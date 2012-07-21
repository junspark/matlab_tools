function d = DistanceToFiber(c, s, quats, sym)
% DISTANCETOFIBER - Distance from points to fiber.
%   VERSION:  $Id: DistanceToFiber.m 136 2009-09-03 13:31:23Z boyce $
%
%   STATUS:  in development
%
%   USAGE:
%
%   d = DistanceToFiber(c, s, quats, sym)
%
%   INPUT:
%
%   c, s are 3 x 1 vectors
%        a crystal and a sample direction, not necessarily normalized;
%        these define the fiber
%   quats is 4 x nq
%            a list of input quaternions
%   sym   is 4 x ns
%            the symmetry group in quaternions
%
%   OUTPUT:
%
%   d is 1 x nq
%        the distance from each quaternion to the (c, s) fiber,
%        accounting for symmetry
%
%   NOTES:
%
%   * Uses functions:  RMatOfQuat, SymHKL, MultMatArray
%   
c = UnitVector(c);
s = UnitVector(s);
%
nquats = size(quats, 2);     % nquats
rmats  = RMatOfQuat(quats);  % 3 x 3 x nquats
csym   = SymHKL(c, sym);     % 3 x m
nhkl   = size(csym, 2);     % m
rc     = MultMatArray(rmats, repmat(csym, [1, 1, nquats])); % 3 x m
                                                            % x nquats
rcdots = s(:)' * reshape(rc, [3, nhkl*nquats]);
rcdots = reshape(rcdots, [nhkl, nquats]);
cosang = max(rcdots, [], 1);
cosang = min(cosang, 1);
d = acos(cosang);
