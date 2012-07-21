function opm = BuildOdfPfMatrix(hkl, mesh, sym, pts, div, invpf, block)
% BuildOdfPfMatrix - Build ODF/PF matrix in pieces.
%   
%   USAGE:
%
%   opm = BuildOdfPfMatrix(hkl, mesh, sym, pts, div, invpf, block)
%
%   INPUT:
%
%   hkl,  mesh, sym, pts, div, invpf 
%         are the same as in `OdfPfMatrix', except that `invpf' is 
%         no longer optional 
%   block is a postive integer, 
%         the number of pole figure points to be processed at a
%         time;  the pole figure points can be processed 
%         independently to control memory allocation
%
%   OUTPUT:
%
%   See `OdfPfMatrix'
%
%   NOTES:
%
%   *  See OdfPfMatrix for further documentation.
%
%   *  This is preferrd to OdfPfMatrix for problems with
%      many divisions per fiber and many pole figure points.
%
npts = size(pts, 2);  
ncol = min(mesh.eqv(1,:) - 1); % assuming eqv nonempty
%
block = min([abs(block) npts]);
%
%-------------------- Compute the matrix.
%
if (block < npts)
  remainder = rem(npts, block);
  nblock    = round((npts - remainder)/block); % ensure integer
else
  remainder = 0;
  nblock    = 1;
end
%
opm = sparse(npts, ncol);
%
for b=1:nblock      % do full blocks
  rmax = b*block;
  rmin = rmax - block + 1; 
  range = rmin:rmax;
  %
  opm(range, :) = OdfPfMatrix(hkl, mesh, sym, pts(:, range), div, invpf);
  %
end
%
if (remainder > 0)
  rmin = npts - remainder + 1;
  range = rmin:npts;
  %
  opm(range, :) = OdfPfMatrix(hkl, mesh, sym, pts(:, range), div, invpf);
  %
end
