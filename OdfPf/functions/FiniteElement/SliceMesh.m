function [sm, smels, smecrd] = SliceMesh(m, p, n)
% SliceMesh - Planar slice through a 3D mesh
%   
%   USAGE:
%
%   sm = SliceMesh(m, p, n)
%
%   INPUT:
%
%   m is a MeshStructure,
%     for a general mesh
%   p is a 3-vector,
%     a point on the slice plane
%   n is a 3-vector,
%     the normal to the slice plane; it does not
%     need to be unit length
%
%   OUTPUT:
%
%   sm is a MeshStructure,
%      on the 2D planar section
%   smels is an integer array,
%      the list of elements containing the nodes in `sm',
%      relative to the original mesh, `m'
%   smecrd is 4 x n,
%      the barycentric (elemental) coordinates of the nodes
%      in `sm' relative to the original mesh, `m'
%
ZEROP = eps;
ZEROM = -ZEROP;
%
p = p(:);
n = UnitVector(n(:));
%
crd = m.crd;
con = m.con;
%
%  First, find elements which intersect.
%
d2p   = DistToPlane(crd, p, n);
edist = d2p(con);
emax  = max(edist);
emin  = min(edist);
els   = find((emax > ZEROM) & (emin< ZEROP));
econ  = con(:, els);
edist = edist(:, els);
%
%  Now find intersection points.
%
ne = length(els);
%
smcrd  = zeros(3, 4*ne); ncrd = 0;
smcon  = zeros(3, 2*ne); nsel = 0;
smels  = zeros(1, 4*ne);
smecrd = zeros(4, 4*ne);
%
for enum = 1:ne
  nodes = econ(:,enum);
  ecrd  = crd(:, nodes);
  ed2p  = edist(:, enum);
  epts  = []; e_els = []; bary = [];
  for p_i = 1:4
    mydi = ed2p(p_i);
    if abs(mydi) < ZEROP
      epts  = [epts, ecrd(:, p_i)];
      e_els = [e_els, els(enum)];
      mybary = [0 0 0 0]'; mybary(p_i) = 1.0;
      bary  = [bary, mybary];
    elseif mydi < ZEROM
      for p_j = 1:4
	mydj = ed2p(p_j);
	if mydj > ZEROP
	  beta  = mydj/(mydj - mydi);
	  newpt = beta*ecrd(:, p_i) + (1 - beta)* ecrd(:, p_j);
	  epts  = [epts, newpt];
	  e_els = [e_els, els(enum)];
	  mybary = [0 0 0 0]';
	  mybary(p_i) = beta;  mybary(p_j) = 1 - beta;
	  bary = [bary, mybary];
	end
      end
    end
    %
    numnew = size(epts, 2);
    if (numnew == 3)
      %
      rng = ncrd+(1:3);
      %
      smcrd(:, rng) = epts;
      smels(rng) = e_els;
      smecrd(:, rng) = bary;
      %
      smcon(:, nsel + 1) = [1; 2; 3] + ncrd;
      %smcrd  = [smcrd, epts];
      %smcon  = [smcon, ([1; 2; 3] + ncrd) ];
      %smels  = [smels, e_els];
      %smecrd = [smecrd, bary];
      ncrd = ncrd + 3;
      nsel = nsel + 1;
      %
    elseif (numnew == 4)
      %
      rng = ncrd+(1:4);
      %
      smcrd(:, rng) = epts;
      smels(rng)     = e_els;
      smecrd(:, rng) = bary;
      %
      smcon(:, nsel + 1) = [1; 2; 3] + ncrd;
      smcon(:, nsel + 2) = [2; 3; 4] + ncrd;
      %smcrd  = [smcrd, epts];
      %smcon  = [smcon, ([1 2; 2 3; 3 4] + ncrd) ];
      %smels  = [smels, e_els];
      %smecrd = [smecrd, bary];
      ncrd = ncrd + 4;
      nsel = nsel + 2;
      %
    end
    %
  end
end
sm = MeshStructure(smcrd(:, 1:ncrd), smcon(:, 1:nsel));
smels = smels(1:ncrd);
smecrd = smecrd(:, 1:ncrd);
%
%--------------------Private Functions----------------------------------
%
function d = DistToPlane(x, p, n)
% DISTTOPLANE - Distance to plane
%   
pr = repmat(p, [1, size(x, 2)]);
d  = n' * (x - pr);
