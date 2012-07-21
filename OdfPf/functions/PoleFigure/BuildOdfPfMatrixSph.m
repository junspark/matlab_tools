function opm = BuildOdfPfMatrixSph(hkl, mesh, sym, pts, div, block, varargin)
% BUILDODFPFMATRIX - Build ODF/PF matrix in pieces.
%   
%   opm = BuildOdfPfMatrix(hkl, mesh, sym, pts, div, invpf, block)
%
%   hkl, mesh, sym, pts, div, invpf are the same as in OdfPfMatrix,
%        except that `invpf' is no longer optional 
%   block is 1 x 1, the number of points to be processed at a time
%
%   Notes:
%
%   *)  This routine is preferred to the basic OdfPfMatrix,
%       which can consume too much memory for a problem with
%       many points per fiber.  
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
    %
    rmax = b*block;
    rmin = rmax - block + 1; 
    range = rmin:rmax;
    %
    fprintf('\n\tBuildOdfPfMatrix: processing points %d:%d', rmin, rmax)
    %
    opm(range, :) = OdfPfMatrixSph(hkl, mesh, sym, pts(:, range), div, varargin{:});
    %
end
%
if (remainder > 0)
    rmin = npts - remainder + 1;
    range = rmin:npts;
    %
    fprintf('\n\tBuildOdfPfMatrix: processing points %d:%d\n', rmin, npts)
    %
    opm(range, :) = OdfPfMatrixSph(hkl, mesh, sym, pts(:, range), div, varargin{:});
    %
end
