function fibs = FindFiber(c, s, ndiv)
% FindFiber - Find fiber for given crystal and sample directions.
%  
%   USAGE:
%
%   fibs = FindFiber(c, s, ndiv)
%
%   INPUT:
%
%   c is 3 x m,
%        an array of crystal directions
%   s is 3 x n,
%        an array of sample directions;
%        if n > 1 then m must be 1, and vice versa
%   ndiv is a postive integer,
%        the number of points to return per fiber
%
%   OUTPUT:
%
%   fibs is 4 x ndiv x nfibs,
%        an array of equally spaced points (quaternions) along each
%        specified fiber; `nfibs' is either `m' or `n' above
%
%   NOTES:
%
%   * The (c, s) fiber is the collection of rotations R
%     such that Rc = s
%
%   * This routine can be called with many crystal directions
%     and a single sample (as for inverse pole figures) or
%     with many sample directions and a single crystal direction
%     (as for pole figures)
%
nc = size(c, 2); % number of crystal directions
ns = size(s, 2); % number of sample directions
%
if (nc > 1) & (ns > 1)
  error('more than one crystal and more than one sample directions')
end
%
npts = max(nc, ns);
%
c = UnitVector(c);
s = UnitVector(s);
%
if (nc > 1)
  h = s;
  s = repmat(s, [1 npts]);
else
  h = c;
  c = repmat(c, [1 npts]);
end
%
cutoff = 1.0e-8;
%
%  We start with a rotation of pi about the midpoint.
%
ax    = c + s;
anrm  = sqrt(dot(ax, ax, 1));
okay  = (anrm > cutoff);
nokay = sum(okay);
%
if (nokay < npts)
  nlsp = null(h');
  hort = nlsp(:, 1);
  ax(:, ~okay) = repmat(hort, [1, npts - nokay]);
end
%
q0 = [zeros(1, npts); UnitVector(ax)];
%
%  Now find rotations which leave h fixed. 
%
phi = (0:(ndiv-1))*((2*pi)/ndiv);
qh  = QuatOfAngleAxis(phi, repmat(h, [1 ndiv]));
%
%  Now multiply.
%
fibs = zeros(4, ndiv, npts);
if (nc > 1)
  for i=1:npts
    fibs(:, :, i) = QuatProd(qh, repmat(q0(:, i), [1 ndiv]));
  end
else
  for i=1:npts
    fibs(:, :, i) = QuatProd(repmat(q0(:, i), [1 ndiv]), qh);
  end
end
