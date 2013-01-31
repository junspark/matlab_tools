function bndry = MeshBoundary(m)
%  MeshBoundary - Generate mesh for boundary. (IN DEVELOPMENT?)
%   
%   USAGE:
%
%   bndry = MeshBoundary(m)
%
%   INPUT:
%
%   m is a MeshStructure
%        The mesh to operate on.
%
%   OUTPUT:
%
%   bndry is a MeshStructure
%         A mesh on the boundary of `m'.  It has additional
%         structure fields of 'parent_nodes', 'parent_elems'
%         and 'isboundary' [=1].
%
szmc = size(m.con); 
%
%  Make list of all element faces, including multiplicity.
%
esurfs = m.etype.surfs;
szes   = size(esurfs);
nallf  = szes(2)*szmc(2);
f      = reshape(m.con(m.etype.surfs, :), [szes(1), nallf]);
%
%  Eliminate multiples.
%
[fsort,    f_ind]      = sort(f);                  [tmp, f_ind_inv]  = sort(f_ind);
[allfaces, fs_ind]     = sortrows(fsort');         [tmp, fs_ind_inv] = sort(fs_ind);
[ufaces,   f_of, a_of] = unique(allfaces, 'rows'); [tmp, f_of_inv]   = sort(f_of);
%
%  Count multiplicity to determine surface elements.
%
nfaces = size(ufaces, 1);
wunz   = ones(1, nallf);
%
mult = sparse(wunz, a_of, wunz, 1, nfaces);
mult = full(mult);   
%
if any(mult > 2)
  error('multiplicity > 3 for element surface')
end
%
%  Find indices in original list of all faces, so that original order is
%  preserved.
%
srfels = find(mult == 1);
%
fnums = fs_ind(f_of(srfels));
bndry_con = f(:, fnums);
%
%  Eliminate unused indices.
%
npnums = unique(bndry_con(:)); numnodes = length(npnums);
newnum = sparse(ones(1, numnodes), npnums, 1:numnodes);
bndry_coords = m.crd(:, npnums);
bndry_con    = full(newnum(bndry_con));
bndry = MeshStructure(bndry_coords, bndry_con, [], ...
		      'ElementType', m.etype.surftype);
%
bndry.isboundary   = 1;
bndry.parent_nodes = npnums(:)';
bndry.parent_elems = 1 + floor((fnums(:)' - 1)./szes(2));
