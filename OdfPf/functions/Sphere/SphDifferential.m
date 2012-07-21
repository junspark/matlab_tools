function diff = SphDifferential(mesh, refpts)
% SphDifferential - Compute differential of mapping to sphere.
%   
%   USAGE:
%
%   diff = SphDifferential(mesh, refpts)
%   
%   INPUT:
%
%   mesh   is a mesh,
%          on a sphere of any dimension
%   refpts is d x n, 
%          a list of points in the reference element,
%          usually the quadrature points, given in barycentric
%          coordinates
%
%   OUTPUT:
%
%   diff is d x (d-1) x nq, 
%        a list of tangent vectors at each reference
%        point for each element; nq is the number of global 
%        quadrature points, that is n x ne, where ne is the
%        number of elements
%        
crd = mesh.crd;
con = mesh.con;
%eqv = mesh.eqv; % equivalence not used here
%
[dr, nr] = size(refpts);
[dc, nn] = size(crd);
[d,  ne] = size(con);
%
if (dc ~= d)
  error('dimension mismatch: coordinates and elements')
end
%
if (dr ~= d)
  error('dimension mismatch: reference points and elements')
end
%
%  First compute tangent vectors at intermediate
%  mapping to inscribed simplex.  Make a copy for 
%  each quadrature point.
%
dm1 = d - 1;
simp_dx = crd(:, con(1:dm1, :)) - ...
	reshape(repmat(crd(:, con(d, :)), [dm1, 1]), ...
		[d, dm1*ne]);
stmp = reshape(simp_dx, [d*dm1 ne]);
stmp = reshape(repmat(stmp, [nr 1]), [d dm1*nr*ne]);
pts  = reshape(SpreadRefPts(mesh, refpts), [d nr*ne]);
ptmp = reshape(repmat(pts, [dm1 1]), [d dm1*nr*ne]);
%
%  Now apply radial projection to sphere.
%
nrmp = 1.0 ./ sqrt(dot(ptmp, ptmp, 1));
fact = repmat((dot(stmp, ptmp, 1).*nrmp.*nrmp), [d 1]);
diff = repmat(nrmp, [d 1]) .* (stmp - (fact.*ptmp));
diff = reshape(diff, [d dm1 nr*ne]);
