function mnew = MeshSort(mold, dir)
% MeshSort - Renumber mesh by spatial sorting.
%   
%   USAGE:
%
%   mnew = MeshSort(mold, dir)
%
%   INPUT:
%
%   mold is a MeshStructure
%
%   dir is an integer (or an integer vector)
%        It gives the spatial direction of the sort.  See `sortrows'
%        for the role of `dir'.
%
%   OUTPUT:
%
%   mnew is a MeshStructure
%        It is the new mesh with coordinates and elements
%        reordered according to `dirs'.
%
%   NOTES:
%
%   * Does not handle equivalences.
%   * Should return indexing arrays for associated data
%
msize = size(mold.con);
mnew  = mold;  % preserve element type
%
%  First, sort nodes and modify connectivity to match.
%
[xsrt, xind] = sortrows(mold.crd', dir);
mnew.crd = xsrt';
[tmp, invord] = sort(xind);
mnew.con = invord(mold.con);
%
%  Second, sort connectivity.
%
tmp = reshape(mnew.crd(dir(1), mnew.con), msize);
tmp = min(tmp);
[tmps, reord] = sort(tmp);
mnew.con = mnew.con(:, reord);
%
%  Remove regular data fields if present, since sorting will
%  alter the order.
%
if isfield(mnew, 'ndiv')
  rmfield(mnew, {'ndiv', 'transform'})
end
%
