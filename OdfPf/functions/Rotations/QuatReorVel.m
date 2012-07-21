function qVel = QuatReorVel(q, w)
% QuatReorVel - Map axial vectors to quaternion tangent vectors
%   
%   VERSION:  $Id$
%
%   STATUS:  in development
%
%   USAGE:
%
%   qVel = QuatReorVel(q, v)
%
%   INPUT:
%
%   q is 4 x n1 x n2 x ...
%        a list of quaternions
%   w is 3 x n1 x n2 x ...
%        a list of axial vectors of skew matrices W, such that
%        R_dot = W R, where w = vec(W)
%
%   OUTPUT:
%
%   qVel is 4 x n1 x n2 x ...
%
%   NOTES:
%
%   * Mapping is:  w |->  (1/2) (-w . qv, qs w + w X qv)
%
sq = size(q); sw = size(w);
qvec = reshape(q(2:4, :), sw);
%
qVel = 0.5 * [-dot(w, qvec); ...
	      reshape(repmat(q(1, :), [3 1]), sw) .* w + ...
	      cross(w, qvec) ];
